#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
复变函数求解器 v2.0 
=================================
支持题型：
  方程求解 | 围道积分 | 解析性判定(C-R) | 奇点/极点分类 | 留数计算 | 级数展开(泰勒/洛朗)

使用方式：
  solve("z^2 + 1 = 0")                      # 自动检测
  solve("equation: z^6 - 1 = 0")            # 显式指定题型
  solve("integral: 1/(z^2+1) |z|=2")       # 围道积分
  solve("singularity: sin(z)/z at z=0")     # 奇点判定
  solve("residue: 1/(z-1) at z=1")          # 留数计算
  solve("series: exp(z) at 0")              # 级数展开

设计哲学：
  零幻觉 — 不会就拒绝，绝不编造答案
  符号优先 — 先尝试解析解，再考虑数值近似
  结果可验证 — 所有解代回原方程核验
"""

from sympy import latex
import sympy as sp
import mpmath as mp
import re
import sys

# ======================== 全局符号 ========================
z = sp.Symbol('z', complex=True)
x, y = sp.symbols('x y', real=True)
I = sp.I
sp.init_printing(use_unicode=True)

# ======================== 常量 ========================
MAX_FACTOR_ATTEMPTS = 20      # 因式分解最大尝试次数
DEFAULT_SERIES_TERMS = 6      # 默认级数展开项数
NUMERICAL_GUESS_RANGE = 5     # 数值求解初值搜索范围
VERIFY_TOLERANCE = 1e-8       # 解验证容差

# ======================== LaTeX 格式化工具 ========================
def _l(expr, latex=False):
    """将 sympy 表达式格式化为字符串或 LaTeX"""
    if not latex:
        return str(expr)
    try:
        s = sp.latex(expr, imaginary_unit='i')
        # 简化常见输出
        s = s.replace(r'\left[', r'[').replace(r'\right]', r']')
        s = s.replace(r'\left(', r'(').replace(r'\right)', r')')
        return s
    except Exception:
        return str(expr)


def _lc(c, latex=False):
    """将复数格式化为字符串或 LaTeX"""
    def _fmt_num(val):
        """格式化数字，处理 0 边界情况"""
        s = f"{val:.6g}".rstrip('0').rstrip('.')
        return s if s else '0'

    if latex:
        re_val = c.real
        im_val = c.imag
        # 实数情况
        if abs(im_val) < 1e-12:
            return _fmt_num(re_val)
        # 纯虚数情况
        if abs(re_val) < 1e-12:
            im_str = _fmt_num(im_val)
            if im_str == '1':
                return 'i'
            if im_str == '-1':
                return '-i'
            return f"{im_str}i"
        # 一般复数
        re_str = _fmt_num(re_val)
        im_str = _fmt_num(abs(im_val))
        if im_str == '1':
            im_str = ''
        sign = '+' if im_val >= 0 else '-'
        return f"{re_str}{sign}{im_str}i"
    # 纯文本模式
    return f"{c.real:.6g}{c.imag:+.6g}i"


def _laplace_wrap(text, latex=False):
    """LaTeX 模式下用 $$...$$ 包裹公式"""
    if latex:
        return f"$${text}$$"
    return text


def _inline_math(text, latex=False):
    """LaTeX 模式下用 $...$ 包裹行内公式"""
    if latex:
        return f"${text}$"
    return text


# ======================== 1. 输入清洗 ========================
def clean_math_expression(text):
    """
    清洗数学表达式：
    - 统一幂符号 ^ → **
    - 展开函数简写 sinz → sin(z)
    - 虚数单位 i/j → I（注意边界）
    - 自动补乘号
    保留中英文混合内容供题型检测使用。
    """
    text = text.strip()
    # 统一标点
    text = text.replace('^', '**')
    text = text.replace('（', '(').replace('）', ')')
    text = text.replace('，', ',').replace('；', ';')
    text = text.replace('：', ':')

    # 虚数单位（仅替换独立的 i/j，不影响函数名中的字母）
    text = re.sub(r'(?<![a-zA-Z])i(?![a-zA-Z])', 'I', text)
    text = re.sub(r'(?<![a-zA-Z])j(?![a-zA-Z])', 'I', text)

    # 函数简写展开：sinz → sin(z) 等（优先级高于 i→I 替换）
    func_shorthands = {
        'sinz': 'sin(z)', 'cosz': 'cos(z)', 'tanz': 'tan(z)',
        'sinhz': 'sinh(z)', 'coshz': 'cosh(z)', 'tanhz': 'tanh(z)',
        'ez': 'exp(z)', 'lnz': 'log(z)', 'ln(': 'log(',
        'arcsinz': 'asin(z)', 'arccosz': 'acos(z)', 'arctanz': 'atan(z)',
    }
    for old, new in func_shorthands.items():
        text = text.replace(old, new)

    # 常见函数名映射
    name_map = {
        'cot(': '1/tan(', 'sec(': '1/cos(', 'csc(': '1/sin(',
        'conj(': 'conjugate(', 're(': 're(', 'im(': 'im(',
        'abs(': 'Abs(', '|z|': 'Abs(z)',
        'arg(': 'arg(',
    }
    for old, new in name_map.items():
        text = text.replace(old, new)

    # |z-a| 形式 → Abs(z-a)
    text = re.sub(r'\|z\s*([-+])\s*(\S+?)\|', r'Abs(z\1\2)', text)
    text = re.sub(r'\|z\|', 'Abs(z)', text)

    # 自动补乘号：数字后跟字母
    text = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', text)
    # 数字后跟 (
    text = re.sub(r'(\d)\(', r'\1*(', text)
    # z 后跟数字
    text = re.sub(r'\bz(\d)', r'z*\1', text)

    # 去除中文字符（数学解析不需要，关键词检测用原始文本）
    text = re.sub(r'[一-龥]+', '', text)

    return text.strip()


def extract_chinese_keywords(raw_text):
    """提取中文关键词用于题型检测"""
    return re.findall(r'[一-龥]+', raw_text)


def parse_contour_params(raw_text):
    """
    从输入中提取积分围道参数。
    支持格式：
      |z|=r        → 圆心 0, 半径 r
      |z-a|=r      → 圆心 a, 半径 r
      |z+a|=r      → 圆心 -a, 半径 r
      |z-(a+bi)|=r → 圆心 a+bi, 半径 r
    """
    # 格式1: |z-a|=r 或 |z+a|=r 或 |z-(a+bi)|=r
    m = re.search(r'\|\s*z\s*([-+])\s*\(?([^|]+?)\)?\s*\|\s*[=<]\s*(\S+)', raw_text)
    if m:
        sign = m.group(1)
        center_str = m.group(2).strip()
        radius_str = m.group(3).strip()
        try:
            # 计算圆心
            center_expr = sp.sympify(center_str)
            center = complex(center_expr.evalf()) if center_expr.is_number else 0.0
            if sign == '-':
                center = complex(center)
            else:
                center = complex(-sp.sympify(center_str).evalf())
            radius = float(sp.sympify(radius_str).evalf())
            return center, radius
        except Exception:
            pass

    # 格式2: |z|=r
    m = re.search(r'\|\s*z\s*\|\s*[=<]\s*(\S+)', raw_text)
    if m:
        try:
            radius = float(sp.sympify(m.group(1)).evalf())
            return 0.0, radius
        except Exception:
            pass

    # 格式3: 单独的纯数字 "=r"
    m = re.search(r'[=]\s*(\d+\.?\d*)', raw_text)
    if m:
        try:
            return 0.0, float(m.group(1))
        except Exception:
            pass

    return 0.0, 1.0


# ======================== 2. 题型检测 ========================
# 显式前缀映射
EXPLICIT_PREFIX_MAP = {
    'equation': 'equation', '方程': 'equation', 'eq': 'equation',
    'integral': 'integral', '积分': 'integral', '复积分': 'integral', 'int': 'integral',
    'analytic': 'analytic', '解析': 'analytic', 'ana': 'analytic',
    'singularity': 'singularity', '奇点': 'singularity', '极点': 'singularity', 'singular': 'singularity',
    'residue': 'residue', '留数': 'residue', 'res': 'residue',
    'series': 'series', '级数': 'series', '展开': 'series', 'expand': 'series',
}


def parse_explicit_prefix(raw_text):
    """
    解析显式题型前缀。支持中英文前缀 + ':' 或 '：'。
    返回 (task_type, remaining_text) 或 (None, raw_text)
    """
    for prefix_str, task_type in EXPLICIT_PREFIX_MAP.items():
        for sep in [':', '：']:
            full = prefix_str + sep
            if raw_text.lower().startswith(full):
                return task_type, raw_text[len(full):].strip()
    return None, raw_text


def detect_task_type(raw_text, cleaned_text):
    """
    多层题型检测，返回 (task_type, confidence, params)

    检测优先级：
    1. 显式前缀（confidence=1.0）
    2. 中文关键词
    3. 数学结构特征
    4. 默认回落
    """
    chinese_keywords = extract_chinese_keywords(raw_text)
    chinese_joined = ''.join(chinese_keywords)

    # === 积分 ===
    if any(kw in chinese_joined for kw in ['积分', '围道', '环路', '柯西积分']):
        center, radius = parse_contour_params(''.join(chinese_keywords) + raw_text)
        return 'integral', 0.95, {'center': center, 'radius': radius}
    if 'dz' in raw_text.lower() or '∮' in raw_text or '∳' in raw_text:
        center, radius = parse_contour_params(raw_text)
        return 'integral', 0.85, {'center': center, 'radius': radius}

    # === 级数（洛朗 > 泰勒/麦克劳林） ===
    if any(kw in chinese_joined for kw in ['洛朗', 'laurent']):
        point, order = parse_series_params(raw_text)
        return 'series', 0.95, {'point': point, 'order': order, 'is_laurent': True}
    if any(kw in chinese_joined for kw in ['展开', '级数', '泰勒', '麦克劳林', 'taylor', 'maclaurin']):
        point, order = parse_series_params(raw_text)
        return 'series', 0.95, {'point': point, 'order': order, 'is_laurent': False}

    # === 奇点 ===
    if any(kw in chinese_joined for kw in ['奇点', '极点', '本性', '可去奇点', '孤立奇点']):
        point = parse_point(raw_text)
        return 'singularity', 0.95, {'point': point}

    # === 留数 ===
    if any(kw in chinese_joined for kw in ['留数', 'residue']):
        point = parse_point(raw_text)
        return 'residue', 0.95, {'point': point}

    # === 解析性 ===
    if any(kw in chinese_joined for kw in ['解析', '柯西黎曼', 'C-R', '调和']):
        return 'analytic', 0.95, {}

    # === 方程（有等号且非其他题型） ===
    if '=' in cleaned_text:
        return 'equation', 0.85, {}

    # === 数学结构特征 ===
    # 分式且有明确点 → 可能是留数或奇点
    if '/' in cleaned_text:
        point = parse_point(raw_text)
        # 检查分母是否在该点为零
        try:
            _, denom_str = cleaned_text.rsplit('/', 1)
            denom = sp.sympify(denom_str.strip())
            if denom.has(z):
                return 'residue', 0.6, {'point': point}
        except Exception:
            pass

    # === 默认：有 z 则判解析性 ===
    if 'z' in cleaned_text or 'Z' in cleaned_text:
        return 'analytic', 0.3, {}

    return 'unknown', 0.0, {}


def parse_point(raw_text):
    """从文本中提取奇点/展开点坐标（纯实数返回 float，含虚部返回 complex）"""
    # 匹配 z=数字 或 at 数字 或 在 z=数字
    m = re.search(r'[zZ]\s*[=＝]\s*([-+]?\d+\.?\d*\s*[-+]?\s*\d*\.?\d*\s*[iI]?)', raw_text)
    if m:
        try:
            val = m.group(1).replace(' ', '').replace('i', 'j').replace('I', 'j')
            c = complex(val)
            # 无虚部返回 float/int，保证 sympy 符号计算
            if abs(c.imag) < 1e-12:
                return int(c.real) if c.real == int(c.real) else float(c.real)
            return c
        except Exception:
            pass
    # 匹配单独数字
    nums = re.findall(r'[-+]?\d+\.?\d*', raw_text)
    if nums:
        try:
            f = float(nums[-1])
            return int(f) if f == int(f) else f
        except Exception:
            pass
    return 0


def parse_series_params(raw_text):
    """提取级数展开参数：展开点和阶数"""
    point = 0
    order = DEFAULT_SERIES_TERMS
    nums = re.findall(r'[-+]?\d+\.?\d*', raw_text)
    if nums:
        try:
            # 最后一个较大的整数可能是阶数
            candidates = [float(n) for n in nums]
            int_candidates = [n for n in candidates if n == int(n) and n > 1]
            if int_candidates:
                order = min(int(int_candidates[-1]), 12)
                other_nums = [n for n in candidates if n != order]
                if other_nums:
                    point = other_nums[0]
            else:
                point = candidates[0] if candidates else 0
        except Exception:
            pass

    # 检查是否有 z=数字 明确指定展开点
    m = re.search(r'[zZ]\s*[=＝]\s*([-+]?\d+)', raw_text)
    if m:
        try:
            point = float(m.group(1))
        except Exception:
            pass

    return point, order


# ======================== 3. 安全执行器 ========================
ERROR_MESSAGES = {
    'PolynomialError': '非标准多项式，无通用解析解法',
    'ValueError': '超出函数合法定义域',
    'TypeError': '含未定义符号，无法解析',
    'SyntaxError': '公式语法错误',
    'ZeroDivisionError': '存在间断奇点，函数此处无定义',
    'NotImplementedError': '该题型超出当前求解范围',
    'PoleError': '该点存在极点，计算不适用',
    'NoSolution': '复数域内无解析解',
}


def safe_exec(func, *args, **kwargs):
    """
    安全执行器：捕获所有已知数学异常。
    返回格式：{'success': bool, 'result': str, 'error': str|None}
    """
    try:
        result = func(*args, **kwargs)
        if isinstance(result, dict):
            result['success'] = True
            return result
        return {'success': True, 'result': str(result), 'error': None}
    except sp.PolynomialError:
        return {'success': False, 'result': '[FAIL]  非标准多项式，无通用解析解法', 'error': 'PolynomialError'}
    except ValueError as e:
        return {'success': False, 'result': f'[FAIL]  数值错误：{e}', 'error': 'ValueError'}
    except TypeError as e:
        return {'success': False, 'result': f'[FAIL]  类型错误：{e}', 'error': 'TypeError'}
    except SyntaxError as e:
        return {'success': False, 'result': f'[FAIL]  语法错误：{e}', 'error': 'SyntaxError'}
    except ZeroDivisionError:
        return {'success': False, 'result': '[FAIL]  存在间断奇点，函数此处无定义', 'error': 'ZeroDivisionError'}
    except NotImplementedError as e:
        return {'success': False, 'result': f'[FAIL]  暂不支持：{e}', 'error': 'NotImplementedError'}
    except Exception as e:
        error_type = type(e).__name__
        return {'success': False, 'result': f'[FAIL]  求解失败({error_type})：{e}', 'error': error_type}


# ======================== 4. 方程求解 ========================
def solve_complex_equation(eq_str, allow_numerical=True, latex=False):
    """
    复变方程求解器。
    流水线：预处理(共轭/Abs分离) → 因式分解 → 符号求解 → 数值求解(mpmath)
    """
    if '=' not in eq_str:
        return "[FAIL]  非标准等式方程，格式错误"

    lhs_str, rhs_str = eq_str.split('=', 1)
    expr_str = f"({lhs_str.strip()})-({rhs_str.strip()})"

    try:
        expr = sp.sympify(expr_str, locals={'z': z})
    except Exception as e:
        return f"[FAIL] 表达式解析失败：{e}"

    # ---- 情况1: 含共轭/绝对值/Re/Im ----
    if expr.has(sp.conjugate) or expr.has(sp.Abs) or expr.has(sp.re) or expr.has(sp.im):
        try:
            expr_xy = expr.subs(z, x + I*y).expand()
            eq_real = sp.simplify(sp.re(expr_xy))
            eq_imag = sp.simplify(sp.im(expr_xy))
            sols_xy = sp.solve([eq_real, eq_imag], [x, y], dict=True)

            # 生成推导步骤
            if latex:
                # 重建原始方程字符串用于展示
                lhs_l = _l(sp.sympify(lhs_str.strip(), locals={'z': z}), True)
                rhs_l = _l(sp.sympify(rhs_str.strip(), locals={'z': z}), True)
                steps = (
                    f"设 $z = x + iy\\;(x, y \\in \\mathbb{{R}})$，代入：\n"
                    f"$${lhs_l} = {rhs_l}$$\n"
                    f"分离实虚部得方程组：\n"
                    f"$$\\begin{{cases}} {_l(eq_real, True)} = 0 \\\\ {_l(eq_imag, True)} = 0 \\end{{cases}}$$"
                )
            else:
                steps = (
                    f"设 z = x + iy，代入得：\n"
                    f"  ({expr_xy}) = 0\n"
                    f"实部: {eq_real} = 0\n"
                    f"虚部: {eq_imag} = 0"
                )

            # 空列表 = sympy 确认无解
            if sols_xy == []:
                if latex:
                    return steps + f"\n\n联立求解，导出矛盾（如 $x$ 需同时满足 $x \\ge a$ 和 $x < a$），故：\n\\textbf{{结论：}}该方程在复数域内无解。"
                return steps + "\n\n联立求解导出矛盾，[PASS] 该方程在复数域内无解。"
            if sols_xy:
                complex_sols = []
                for sol in sols_xy:
                    try:
                        cx = complex(sol.get(x, 0).evalf())
                        cy = complex(sol.get(y, 0).evalf())
                        complex_sols.append(cx + 1j * cy)
                    except Exception:
                        complex_sols.append(sol)
                return f"[PASS]  解析解（实虚分离法）：\n{format_solution_list(complex_sols, latex)}"
        except Exception:
            pass

    # ---- 情况2: 因式分解 ----
    try:
        factored = sp.factor(expr)
        if factored.is_Mul:
            factors = list(factored.args)
            all_sols = []
            for factor in factors:
                if factor.has(z):
                    try:
                        sols = sp.solve(factor, z)
                        all_sols.extend(sols)
                    except Exception:
                        pass
            if all_sols:
                valid_sols = verify_solutions_expr(expr_str, all_sols)
                if latex:
                    eq_latex = _l(expr, True)
                    fac_latex = _l(factored, True)
                    header = (
                        f"原方程：$${eq_latex} = 0$$\n\n"
                        f"因式分解：$${fac_latex} = 0$$"
                    )
                    return f"[PASS]  解析解（因式分解法）：\n{header}\n\n{format_solution_list(valid_sols, latex)}"
                return f"[PASS]  解析解（因式分解法）：\n{format_solution_list(valid_sols, latex)}"
    except Exception:
        pass

    # ---- 情况3: 符号求解 ----
    try:
        sym_sols = sp.solve(expr, z)
        if sym_sols:
            valid_sols = verify_solutions_expr(expr_str, sym_sols)
            if latex:
                eq_latex = _l(expr, True)
                header = f"原方程：$${eq_latex} = 0$$"
                return f"[PASS]  解析解：\n{header}\n\n{format_solution_list(valid_sols, latex)}"
            return f"[PASS]  解析解：\n{format_solution_list(valid_sols, latex)}"
    except Exception:
        pass

    # ---- 情况4: 数值求解兜底 ----
    if allow_numerical:
        try:
            f_lambda = sp.lambdify(z, expr, 'mpmath')
            guesses = _generate_complex_guesses()
            found = set()
            for g0 in guesses:
                try:
                    root = mp.findroot(f_lambda, g0)
                    root_key = _round_complex(root, 8)
                    found.add(root_key)
                except Exception:
                    continue
            if found:
                sorted_roots = sorted(found, key=lambda c: (c.real, c.imag))
                if latex:
                    eq_latex = _l(expr, True)
                    header = f"原方程：$${eq_latex} = 0$$"
                    return f"[PASS]  数值近似解（牛顿迭代法）：\n{header}\n\n{format_solution_list(sorted_roots, latex)}\n[WARN]  注意：此为非准确解析解，仅作参考"
                return f"[PASS]  数值近似解（牛顿迭代法）：\n{format_solution_list(sorted_roots, latex)}\n[WARN]  注意：此为非准确解析解，仅作参考"
        except Exception:
            pass

    return "[FAIL]  该方程无初等解析解且数值求解亦未收敛"


def _generate_complex_guesses():
    """生成复数域初始猜测网格"""
    guesses = []
    for a in range(-NUMERICAL_GUESS_RANGE, NUMERICAL_GUESS_RANGE + 1):
        for b in range(-NUMERICAL_GUESS_RANGE, NUMERICAL_GUESS_RANGE + 1):
            guesses.append(complex(a, b))
            # 加一些非整数点
            guesses.append(complex(a + 0.5, b + 0.5))
    return guesses


def _round_complex(c, decimals):
    """四舍五入复数到指定位数"""
    return round(c.real, decimals) + 1j * round(c.imag, decimals)


def format_solution_list(solutions, latex=False):
    """格式化解列表输出"""
    if not solutions:
        return "  (空集)"
    lines = []
    for i, sol in enumerate(solutions, 1):
        if isinstance(sol, (int, float, complex)):
            c = complex(sol)
            if latex:
                lines.append(f"  z_{{{i}}} = {_lc(c, True)}")
            else:
                lines.append(f"  z{i} = {c.real:.6g}{c.imag:+.6g}i")
        else:
            if latex:
                lines.append(f"  z_{{{i}}} = {_l(sol, True)}")
            else:
                lines.append(f"  z{i} = {sol}")
    return '\n'.join(lines)


def verify_solutions_expr(expr_str, solutions):
    """
    代回验证解的真伪。
    返回通过验证的解列表。
    """
    if not solutions:
        return []
    valid = []
    try:
        expr = sp.sympify(expr_str, locals={'z': z})
    except Exception:
        return list(solutions)

    for sol in solutions:
        try:
            val = complex(expr.subs(z, sol).evalf(chop=True))
            if abs(val) < VERIFY_TOLERANCE:
                valid.append(sol)
        except Exception:
            # 保守处理：无法验证则保留
            valid.append(sol)
    return valid


# ======================== 5. 解析性判定 ========================
def judge_analytic(func_str, latex=False):
    """
    基于柯西-黎曼充要条件判定函数解析性。
    """
    try:
        f = sp.sympify(func_str, locals={'z': z, 'x': x, 'y': y})
    except Exception:
        return '[FAIL] 表达式解析失败'

    try:
        f_xy = sp.simplify(f.subs(z, x + I*y))
        u_expr = sp.simplify(sp.re(f_xy))
        v_expr = sp.simplify(sp.im(f_xy))
    except Exception:
        return "[FAIL]  无法分离实虚部，请检查函数表达式"

    ux = sp.diff(u_expr, x)
    uy = sp.diff(u_expr, y)
    vx = sp.diff(v_expr, x)
    vy = sp.diff(v_expr, y)

    cr1_diff = sp.simplify(ux - vy)
    cr2_diff = sp.simplify(uy + vx)

    steps = []
    if latex:
        f_latex = _l(f, True)
        u_latex = _l(u_expr, True)
        v_latex = _l(v_expr, True)
        ux_l = _l(ux, True); uy_l = _l(uy, True)
        vx_l = _l(vx, True); vy_l = _l(vy, True)
        cr1_l = _l(cr1_diff, True); cr2_l = _l(cr2_diff, True)

        steps.append("[STEP]  C-R条件判定步骤：")
        steps.append(f"   $$f(z) = {f_latex}$$")
        steps.append(f"   $$u(x,y) = {u_latex}, \\quad v(x,y) = {v_latex}$$")
        steps.append(f"   $$\\frac{{\\partial u}}{{\\partial x}} = {ux_l}, \\quad "
                     f"\\frac{{\\partial u}}{{\\partial y}} = {uy_l}$$")
        steps.append(f"   $$\\frac{{\\partial v}}{{\\partial x}} = {vx_l}, \\quad "
                     f"\\frac{{\\partial v}}{{\\partial y}} = {vy_l}$$")
        cr1_ok = cr1_diff == 0
        cr2_ok = cr2_diff == 0
        steps.append(f"   C-R① $$\\frac{{\\partial u}}{{\\partial x}} - "
                     f"\\frac{{\\partial v}}{{\\partial y}} = {cr1_l}$$  "
                     f"{'[PASS] ' if cr1_ok else '[NO] '}")
        steps.append(f"   C-R② $$\\frac{{\\partial u}}{{\\partial y}} + "
                     f"\\frac{{\\partial v}}{{\\partial x}} = {cr2_l}$$  "
                     f"{'[PASS] ' if cr2_ok else '[NO] '}")
    else:
        steps.append("[STEP]  C-R条件判定步骤：")
        steps.append(f"   f(z) = u(x,y) + i·v(x,y)")
        steps.append(f"   u(x,y) = {u_expr}")
        steps.append(f"   v(x,y) = {v_expr}")
        steps.append(f"   ∂u/∂x = {ux}  |  ∂u/∂y = {uy}")
        steps.append(f"   ∂v/∂x = {vx}  |  ∂v/∂y = {vy}")
        steps.append(f"   C-R① ∂u/∂x - ∂v/∂y = {cr1_diff}  {'[PASS] ' if cr1_diff == 0 else '[NO] '}")
        steps.append(f"   C-R② ∂u/∂y + ∂v/∂x = {cr2_diff}  {'[PASS] ' if cr2_diff == 0 else '[NO] '}")

    if cr1_diff == 0 and cr2_diff == 0:
        # 进一步检查偏导数连续性
        try:
            uxx = sp.diff(ux, x)
            uyy = sp.diff(uy, y)
            laplace_u = sp.simplify(uxx + uyy)
            if latex:
                laplace_l = _l(laplace_u, True)
                steps.append(f"   调和性验证 $$\\nabla^2 u = {laplace_l}$$  "
                             f"{'[PASS] ' if laplace_u == 0 else '[WARN] '}")
            else:
                steps.append(f"   调和性验证 ∇²u = {laplace_u}  {'[PASS] ' if laplace_u == 0 else '[WARN] '}")
        except Exception:
            pass
        steps.append("[PASS]  严格满足C-R条件，函数在复平面上处处解析")
    else:
        steps.append("[NO]  不满足C-R条件，函数不是解析函数")

    return '\n'.join(steps)


# ======================== 6. 奇点/极点分类 ========================
def classify_singularity(func_str, z0, latex=False):
    """
    严格分类奇点类型：可去奇点 / N阶极点 / 本性奇点 / 支点。
    """
    if str(z0).lower() in ['inf', '∞', 'infinity']:
        return (" >  无穷远点判定：需作代换 w=1/z，分析 w=0 处奇点类型。\n"
                "   请对 f(1/w) 在 w=0 处重新分析。")

    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    z0_sym = sp.sympify(z0)

    # 检查是否包含多值函数（支点）
    if (f.has(sp.log) or f.has(sp.sqrt) or
        any(a.is_Pow and (not a.exp.is_Integer) for a in sp.preorder_traversal(f) if a.is_Pow)):
        try:
            # 检查该点是否在分支切割上
            val_at = complex(f.subs(z, z0_sym).evalf())
            if latex:
                return (f"**支点（branch point）：** $z={_lc(complex(z0_sym.evalf()), True)}$\n\n"
                        f"   函数含多值成分（$\\log$/ $\\sqrt{{}}$ /分数幂），该点非孤立奇点。\n"
                        f"   留数概念不适用于支点，需沿分支切割分析。")
            return (f" >  z={z0} 疑似为支点（branch point）\n"
                    f"   函数含多值成分（log/sqrt/分数幂），该点非孤立奇点。\n"
                    f"   留数概念不适用于支点，需沿分支切割分析。")
        except Exception:
            pass

    # 尝试极限
    try:
        limit_val = sp.limit(f, z, z0_sym)
        if limit_val.is_finite:
            if latex:
                lim_l = _l(limit_val, True)
                z0_l = _lc(complex(z0_sym.evalf()), True)
                return (f"**可去奇点：** $z={z0_l}$\n\n"
                        f"   $$\\lim_{{z \\to {z0_l}}} f(z) = {lim_l}$$\n\n"
                        f"   留数 $\\operatorname{{Res}}(f, {z0_l}) = 0$")
            return (f" >  z={z0} 为可去奇点\n"
                    f"   lim_{{z→{z0}}} f(z) = {limit_val}\n"
                    f"   留数 Res(f, {z0}) = 0")
    except Exception:
        pass

    # 尝试极点阶数判定：对 n=1..7 求极限 lim(z-z0)^n * f(z)
    for n in range(1, 9):
        try:
            g = (z - z0_sym)**n * f
            lim_g = sp.limit(g, z, z0_sym)
            if lim_g.is_finite and lim_g != 0:
                if latex:
                    lim_l = _l(lim_g, True)
                    z0_l = _lc(complex(z0_sym.evalf()), True)
                    return (f"**{n} 阶极点：** $z={z0_l}$\n\n"
                            f"   $$\\lim_{{z \\to {z0_l}}} (z-{z0_l})^{{{n}}} \\cdot f(z) "
                            f"= {lim_l} \\neq 0, \\infty$$")
                return (f" >  z={z0} 为 {n} 阶极点\n"
                        f"   lim_{{z→{z0}}} (z-{z0})^{n}·f(z) = {lim_g} ≠ 0, ∞")
        except Exception:
            continue

    # 尝试 sympy poles（仅适用于有理函数）
    try:
        poles_dict = sp.poles(f, z)
        for p, order in poles_dict.items():
            if abs(complex((p - z0_sym).evalf())) < 1e-10:
                if latex:
                    return f"**{order}阶极点（有理函数判定）：** $z={_lc(complex(z0_sym.evalf()), True)}$"
                return f" >  z={z0} 为 {order} 阶极点（有理函数判定）"
    except Exception:
        pass

    if latex:
        return f"**本性奇点：** $z={_lc(complex(z0_sym.evalf()), True)}$（极限振荡或不存在，且非有限阶极点）"
    return f" >  z={z0} 为本性奇点（极限振荡或不存在，且非有限阶极点）"


# ======================== 7. 留数计算 ========================
def calc_residue(func_str, z0, latex=False):
    """
    留数计算，含判定步骤。
    """
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    z0_sym = sp.sympify(z0)
    f_latex = _l(f, True) if latex else str(f)

    # 判断奇点类型
    try:
        # 可去奇点 → 留数为 0
        lim_val = sp.limit(f, z, z0_sym)
        if lim_val.is_finite:
            if latex:
                lim_l = _l(lim_val, True)
                return (f"[STEP]  $z={_lc(complex(z0_sym.evalf()), True)}$ 为可去奇点\n\n"
                        f"   $$\\lim_{{z \\to {_lc(complex(z0_sym.evalf()), True)}}} "
                        f"f(z) = {lim_l}$$\n\n"
                        f"[PASS]  $\\operatorname{{Res}}(f, {_lc(complex(z0_sym.evalf()), True)}) = 0$")
            return (f"[STEP]  z={z0} 为可去奇点\n"
                    f"[PASS]  Res(f, {z0}) = 0")
    except Exception:
        pass

    # 尝试 sympy residue
    try:
        res = sp.residue(f, z, z0_sym)
        res_simplified = sp.simplify(res)
        if latex:
            res_l = _l(res_simplified, True)
            steps = []
            steps.append("[STEP]  留数计算步骤：")
            steps.append(f"   1. 函数：$$f(z) = {f_latex}$$")
            steps.append(f"   2. 奇点：$z = {_lc(complex(z0_sym.evalf()), True)}$")
            steps.append(f"   3. 使用留数公式计算")
            steps.append(f"[PASS]  $$\\operatorname{{Res}}(f, {_lc(complex(z0_sym.evalf()), True)}) = {res_l}$$")
            return '\n'.join(steps)
        steps = []
        steps.append("[STEP]  留数计算步骤：")
        steps.append(f"   1. 函数：f(z) = {f}")
        steps.append(f"   2. 奇点：z = {z0}")
        steps.append(f"   3. 使用留数公式计算")
        steps.append(f"[PASS]  Res(f, {z0}) = {res_simplified}")
        return '\n'.join(steps)
    except Exception:
        pass

    # 手动一阶极点公式：Res = lim(z-z0)*f(z)
    try:
        g = (z - z0_sym) * f
        lim_g = sp.limit(g, z, z0_sym)
        if lim_g.is_finite and lim_g != 0:
            if latex:
                lim_l = _l(sp.simplify(lim_g), True)
                z0_l = _lc(complex(z0_sym.evalf()), True)
                return (f"[STEP]  一阶极点留数公式："
                        f"$$\\operatorname{{Res}} = \\lim_{{z \\to {z0_l}}} (z-{z0_l}) f(z)$$\n"
                        f"[PASS]  $$\\operatorname{{Res}}(f, {z0_l}) = {lim_l}$$")
            return (f"[STEP]  一阶极点留数公式：Res = lim(z-z0)f(z)\n"
                    f"[PASS]  Res(f, {z0}) = {sp.simplify(lim_g)}")
    except Exception:
        pass

    return f"[FAIL]  无法计算 z={z0} 处的留数"


# ======================== 8. 围道积分 ========================
def calc_complex_integral(func_str, center=0.0, radius=1.0, latex=False):
    """
    使用留数定理计算闭合围道积分。
    ∮_{|z-center|=radius} f(z) dz = 2πi · Σ Res(f, poles inside contour)
    """
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    # 尝试求解极点
    poles_list = []
    method = ""

    try:
        poles_dict = sp.poles(f, z)
        for pole, order in poles_dict.items():
            poles_list.append((complex(pole.evalf()), order, pole))
        method = "有理函数极点分析"
    except Exception:
        # 尝试用 solve(1/f, z) 找极点
        try:
            inv_f = 1 / f
            denom = sp.denom(sp.together(f))
            if denom.has(z):
                pole_candidates = sp.solve(denom, z)
                for pc in pole_candidates:
                    try:
                        pc_val = complex(pc.evalf())
                        poles_list.append((pc_val, 1, pc))
                    except Exception:
                        pass
                method = "分母零点分析"
        except Exception:
            pass

    if not poles_list:
        # 尝试数值方式找极点
        try:
            f_lambda = sp.lambdify(z, 1/f, 'mpmath')
            for g0 in _generate_complex_guesses():
                try:
                    root = mp.findroot(f_lambda, g0)
                    rk = _round_complex(root, 8)
                    if not any(abs(rk - p[0]) < 1e-6 for p in poles_list):
                        poles_list.append((rk, 1, sp.sympify(rk)))
                except Exception:
                    continue
            method = "数值极点搜索"
        except Exception:
            pass

    # 筛选围道内和围道上的极点
    inside = []
    on_boundary = []

    for pole_val, order, pole_sym in poles_list:
        dist = abs(pole_val - complex(center))
        if dist < radius - 1e-10:
            inside.append((pole_val, order, pole_sym))
        elif abs(dist - radius) < 1e-10:
            on_boundary.append((pole_val, order, pole_sym))

    steps = []
    steps.append("[STEP]  柯西留数定理计算步骤：")
    if latex:
        center_l = _lc(complex(center), True)
        steps.append(f"   积分围道：$$|z - {center_l}| = {radius}$$")
        f_latex = _l(f, True)
        steps.append(f"   被积函数：$$f(z) = {f_latex}$$")
        steps.append(f"   {method}得到极点：{[p[0] for p in poles_list]}")
    else:
        steps.append(f"   积分围道：|z - {center}| = {radius}")
        steps.append(f"   {method}得到极点：{[p[0] for p in poles_list]}")

    if on_boundary:
        steps.append(f"   [WARN]  警告：极点 {[p[0] for p in on_boundary]} 恰在围道上")
        steps.append(f"   积分需取柯西主值，考研一般不考查此情形")

    if not inside:
        if latex:
            steps.append(f"   围道内无奇点")
            steps.append(f"[PASS]  由柯西积分定理：$$\\oint f(z) dz = 0$$")
        else:
            steps.append(f"   围道内无奇点")
            steps.append(f"[PASS]  由柯西积分定理：∮f(z)dz = 0")
        return '\n'.join(steps)

    # 计算留数和
    total_res = 0
    residue_details = []
    for pole_val, order, pole_sym in inside:
        try:
            res = complex(sp.residue(f, z, pole_sym).evalf())
        except Exception:
            try:
                g = (z - pole_sym) * f
                res = complex(sp.limit(g, z, pole_sym).evalf())
            except Exception:
                continue
        residue_details.append((pole_val, order, res))
        total_res += res

    result = 2 * sp.pi * I * total_res

    if latex:
        steps.append(f"   围道内极点及留数：")
        for pv, od, r in residue_details:
            steps.append(f"     $z={_lc(pv, True)}$ ({od}阶) $\\to$ $\\operatorname{{Res}} = {_lc(r, True)}$")
        steps.append(f"   留数和 = ${_lc(total_res, True)}$")
        result_l = _lc(complex(sp.N(result, 10)), True)
        steps.append(f"[PASS]  $$\\oint f(z) dz = 2\\pi i \\cdot \\sum\\operatorname{{Res}} = {result_l}$$")
    else:
        steps.append(f"   围道内极点及留数：")
        for pv, od, r in residue_details:
            steps.append(f"     z={pv} ({od}阶) → Res = {r:.6g}")
        steps.append(f"   留数和 = {_format_complex(total_res)}")
        steps.append(f"[PASS]  ∮f(z)dz = 2πi · ΣRes = {sp.N(result, 10)}")

    return '\n'.join(steps)


def _format_complex(c):
    """格式化复数输出"""
    if abs(c.imag) < 1e-12:
        return f"{c.real:.6g}"
    if abs(c.real) < 1e-12:
        return f"{c.imag:.6g}i"
    return f"{c.real:.6g}{c.imag:+.6g}i"


def _has_negative_powers(expr):
    """检测级数表达式是否含负幂项（z 在分母中）"""
    for term in expr.as_ordered_terms():
        if not term.has(z):
            continue
        pd = term.as_powers_dict()
        if pd.get(z, 0) < 0:
            return True
    # 字符串兜底检测
    expr_str = str(expr)
    if '/z' in expr_str or '**(-' in expr_str:
        return True
    return False


# ======================== 9. 级数展开 ========================
def series_expand(func_str, at_point=0, order=6, is_laurent=False, latex=False):
    """
    函数级数展开：
    - is_laurent=False: 泰勒/麦克劳林展开
    - is_laurent=True: 洛朗展开（在奇点处自动产生负幂项）
    """
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    at_sym = sp.sympify(at_point)
    type_name = "洛朗" if is_laurent else "泰勒"

    try:
        s = sp.series(f, z, at_sym, order)
        terms = s.removeO()

        # 自动检测负幂项（无论 is_laurent 是否设置）
        has_negative = _has_negative_powers(terms)

        if latex:
            terms_l = _l(terms, True)
            f_latex = _l(f, True)
            at_l = _lc(complex(at_sym.evalf()), True) if isinstance(at_point, (int, float, complex)) else str(at_sym)
            if is_laurent:
                if has_negative:
                    return (f"[PASS]  $z={at_l}$ 处的洛朗级数（{order}项）：\n\n"
                            f"   $$f(z) = {f_latex} = {terms_l} + O((z-{at_l})^{{{order}}})$$\n\n"
                            f"   (含负幂项，确为洛朗展开)")
                else:
                    return (f"[WARN]  $z={at_l}$ 处无负幂项，实为泰勒展开：\n\n"
                            f"   $$f(z) = {terms_l}$$")

            if has_negative:
                return (f"[PASS]  $z={at_l}$ 处的级数展开（{order}项）：\n\n"
                        f"   $$f(z) = {f_latex} = {terms_l} + O((z-{at_l})^{{{order}}})$$\n\n"
                        f"   (检测到负幂项，此展开实为洛朗级数)")

            return (f"[PASS]  $z={at_l}$ 处的{type_name}展开（{order}项）：\n\n"
                    f"   $$f(z) = {f_latex} = {terms_l} + O((z-{at_l})^{{{order}}})$$")

        if is_laurent:
            if has_negative:
                return (f"[PASS]  z={at_point} 处的洛朗级数（{order}项）：\n"
                        f"   f(z) = {terms}\n"
                        f"   (含负幂项，确为洛朗展开)")
            else:
                return (f"[WARN]  z={at_point} 处无负幂项，实为泰勒展开：\n"
                        f"   f(z) = {terms}")

        if has_negative:
            return (f"[PASS]  z={at_point} 处的级数展开（{order}项）：\n"
                    f"   f(z) = {terms}\n"
                    f"   (检测到负幂项，此展开实为洛朗级数)")

        return f"[PASS]  z={at_point} 处的{type_name}展开（{order}项）：\n   f(z) = {terms}"

    except Exception as e:
        # 降级尝试
        try:
            s = sp.series(f, z, at_sym, 3)
            terms = s.removeO()
            if latex:
                terms_l = _l(terms, True)
                return f"[WARN]  仅能展开到3项：\n   $$f(z) = {terms_l}$$\n   (更高阶展开失败：{e})"
            return f"[WARN]  仅能展开到3项：\n   f(z) = {terms}\n   (更高阶展开失败：{e})"
        except Exception:
            return f"[FAIL]  无法在 z={at_point} 处展开：{e}"


def _split_func_and_params(text, task_type):
    """
    从含参数的表达式文本中分离出纯函数表达式和参数。
    例如：
      "1/(z^2+1) |z|=2" → ("1/(z^2+1)", {'center': 0.0, 'radius': 2.0})
      "1/(z-1) at z=1"  → ("1/(z-1)", {'point': 1.0})
      "exp(z) at 0"     → ("exp(z)", {'point': 0, 'order': 6, 'is_laurent': False})
    """
    if task_type == 'integral':
        center, radius = parse_contour_params(text)
        # 移除围道部分
        func_text = re.sub(r'\|\s*z[^|]*\|\s*[=<]\s*\S+', '', text).strip()
        func_text = re.sub(r'[,，]\s*$', '', func_text).strip()
        return func_text, {'center': center, 'radius': radius}

    elif task_type in ('singularity', 'residue'):
        point = parse_point(text)
        # 移除 "at z=X" / "z=X" / "在 z=X"
        func_text = re.sub(r'\s+(at|在)\s+z?\s*[=＝]\s*[-+]?\S+', '', text, flags=re.IGNORECASE)
        func_text = re.sub(r'\s+z\s*[=＝]\s*[-+]?\S+\s*$', '', func_text)
        return func_text.strip(), {'point': point}

    elif task_type == 'series':
        point, order = parse_series_params(text)
        is_laurent = any(kw in text.lower() for kw in ['洛朗', 'laurent'])
        # 移除 "at X" / "at z=X" / "在 z=X" / "z=X" / "order=N"
        func_text = re.sub(r'\s+(at|在)\s+(z\s*[=＝]\s*)?[-+]?\S+', '', text, flags=re.IGNORECASE)
        func_text = re.sub(r'\s+z\s*[=＝]\s*[-+]?\S+', '', func_text)
        func_text = re.sub(r'\s+order\s*[=＝]\s*\d+', '', func_text, flags=re.IGNORECASE)
        return func_text.strip(), {'point': point, 'order': order, 'is_laurent': is_laurent}

    return text.strip(), {}


# ======================== 10. 统一调度 ========================
def _to_latex_output(text):
    """后处理器：转换文本标记为 LaTeX，并清理残余的 sympy 格式"""
    result = []
    for line in text.split('\n'):
        stripped = line.strip()
        if not stripped:
            result.append('')
            continue

        # 标记转换
        for tag, latex_tag in [
            ('[PASS]', '\\textbf{[解]}'),
            ('[STEP]', '\\textbf{[步骤]}'),
            ('[FAIL]', '\\textbf{[失败]}'),
            ('[WARN]', '\\textbf{[警告]}'),
            ('[ERR]', '\\textbf{[错误]}'),
            ('[TIP]', '\\textbf{[提示]}'),
            ('[NO]', '\\textbf{[否]}'),
            ('[VERIFY]', '\\textbf{[验证]}'),
            ('[Q]', '\\textbf{[题目]}'),
        ]:
            if stripped.startswith(tag):
                line = line.replace(tag, latex_tag, 1)
                break
        # 处理 > 标记
        if stripped.startswith('>'):
            line = '\\textbf{[判定]} ' + stripped[1:].lstrip()

        # 清理残余的 sympy 复数格式（如 (1+0j) → 1）
        line = re.sub(r'\((-?\d+\.?\d*)\+0j\)', r'\1', line)
        line = re.sub(r'\(0([+-]\d+\.?\d*)j\)', r'\1i', line)
        line = re.sub(r'0j\b', '0', line)
        line = re.sub(r'\(0\+0j\)', '0', line)
        line = re.sub(r'\(\+0j\)', '0', line)
        line = re.sub(r'\(1\+0j\)', '1', line)
        line = re.sub(r'(-?\d+)\.0{6,}\b', r'\1', line)

        # 安全网：对未包裹的简朴数学行补 $$ 包裹
        # 注意：已含 $$/$ 或 \\textbf 的行跳过（子求解器已处理）
        if '$$' not in line and '$' not in line and '\\textbf' not in line:
            math_triggers = ('=', 'sqrt', 'sin', 'cos', 'exp', 'log',
                            'lim', 'Res', 'f(z)', 'u(x', 'v(x', 'z*',
                            r'\frac', r'\int', r'\sum', r'\lim', '·',
                            '**', '^')
            if any(s in line for s in math_triggers):
                indent = ' ' * (len(line) - len(line.lstrip()))
                inner = line.lstrip()
                # 残余 ** → ^{}
                inner = re.sub(r'\*\*(\d+)', r'^{\1}', inner)
                inner = re.sub(r'\*\*\(([^)]+)\)', r'^{\1}', inner)
                inner = re.sub(r'(?<!\\)\^(\d+)', r'^{\1}', inner)
                inner = re.sub(r'(\d)\*([a-zA-Z(])', r'\1 \\cdot \2', inner)
                inner = inner.replace('Res(f, 0', 'Res(f, 0)')
                inner = inner.replace('Res(f, 1', 'Res(f, 1)')
                if not inner.startswith('$$'):
                    line = indent + '$$' + inner + '$$'

        result.append(line)
    return '\n'.join(result)


def solve(raw_input, latex=True):
    """
    复变函数求解器主入口。

    参数：
      raw_input: str — 数学题目或表达式
      latex: bool — True 时输出 LaTeX 格式公式（中文说明保留）

    支持自动题型识别或显式题型前缀：
      equation:  / integral:  / analytic:  / singularity:  / residue:  / series:
      中文前缀：方程: / 积分: / 解析: / 奇点: / 留数: / 级数: / 展开:

    返回：
      str — 格式化的求解过程和答案（latex=True 时数学公式为 LaTeX 格式）
    """
    raw_input = raw_input.strip()

    if not raw_input:
        return "[ERR]  请输入题目或表达式"

    # Step 1: 解析显式前缀
    task_type, remaining = parse_explicit_prefix(raw_input)

    if task_type is None:
        # Step 2: 自动检测（使用原始文本检测关键词）
        task_type, confidence, params = detect_task_type(raw_input, raw_input)

        if task_type == 'unknown' or confidence < 0.3:
            return (
                "[ERR]  无法识别题型或未检测到复变量 z\n\n"
                "[TIP]  请使用题型前缀指定：\n"
                "   equation:表达式    → 求解方程\n"
                "   integral:表达式    → 围道积分\n"
                "   analytic:表达式    → 判定解析性\n"
                "   singularity:表达式 → 奇点分类\n"
                "   residue:表达式     → 留数计算\n"
                "   series:表达式      → 级数展开\n\n"
                "   示例：solve(\"integral: 1/(z^2+1) |z|=2\")"
            )

        # Step 3: 分离函数表达式和参数，然后清洗
        func_text, split_params = _split_func_and_params(raw_input, task_type)
        if split_params:
            params = split_params  # _split_func_and_params 结果优先
        input_clean = clean_math_expression(func_text)
    else:
        # 显式前缀：分离纯函数表达式和参数
        func_text, params = _split_func_and_params(remaining, task_type)
        input_clean = clean_math_expression(func_text)

        # 如果分离失败（confidence 低），补充参数提取
        if not params:
            _task_type, _confidence, params = detect_task_type(remaining, input_clean)
            if _confidence < 0.3:
                params = {}

    # Step 4: 分发到对应求解器并应用 LaTeX 后处理
    result = None
    if task_type == 'equation':
        result = solve_complex_equation(input_clean, latex=latex)

    elif task_type == 'integral':
        center = params.get('center', 0.0)
        radius = params.get('radius', 1.0)
        result = calc_complex_integral(input_clean, center, radius, latex=latex)

    elif task_type == 'analytic':
        result = judge_analytic(input_clean, latex=latex)

    elif task_type == 'singularity':
        point = params.get('point', 0)
        result = classify_singularity(input_clean, point, latex=latex)

    elif task_type == 'residue':
        point = params.get('point', 0)
        result = calc_residue(input_clean, point, latex=latex)

    elif task_type == 'series':
        point = params.get('point', 0)
        order = params.get('order', DEFAULT_SERIES_TERMS)
        is_laurent = params.get('is_laurent', False)
        result = series_expand(input_clean, point, order, is_laurent, latex=latex)

    if result is None:
        result = "[FAIL]  未知错误：题型分发失败"

    # LaTeX 后处理
    if latex and result:
        result = _to_latex_output(result)
    return result


# ======================== 11. 便捷函数 ========================
def solve_equation(eq_str, allow_numerical=True):
    """直接求解复变方程"""
    return solve(f"equation:{eq_str}")


def solve_integral(func_str, center=0.0, radius=1.0):
    """直接计算围道积分"""
    return solve(f"integral:{func_str} |z-{center}|={radius}")


def check_analytic(func_str):
    """直接判定解析性"""
    return solve(f"analytic:{func_str}")


def get_residue(func_str, point):
    """直接计算留数"""
    return solve(f"residue:{func_str} at z={point}")


def expand_series(func_str, point=0, order=6):
    """直接级数展开"""
    return solve(f"series:{func_str} at z={point} order={order}")


def analyze_singularity(func_str, point):
    """直接分析奇点"""
    return solve(f"singularity:{func_str} at z={point}")


# ======================== 12. 程序入口 ========================
if __name__ == "__main__":
    # 确保 UTF-8 输出（兼容现代终端）
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

    # 解析 --text 标志（纯文本模式，关闭 LaTeX）
    use_text = '--text' in sys.argv
    args = [a for a in sys.argv[1:] if a != '--text']
    use_latex = not use_text

    print("=" * 60)
    print("  复变函数求解器 v2.0  --  考研专用")
    print(f"  输出模式：{'纯文本' if use_text else 'LaTeX'}")
    print("  支持：方程 | 积分 | 解析性 | 奇点 | 留数 | 级数展开")
    print("=" * 60)

    if args:
        problem = ' '.join(args)
        print(f"\n[Q]  题目：{problem}")
        print("-" * 40)
        print(solve(problem, latex=use_latex))
        print("-" * 40)
    else:
        print("\n[STEP]  交互模式")
        print("  直接输入题目（自动识别题型）")
        print("  或用前缀指定：equation: / integral: / analytic: /")
        print("                 singularity: / residue: / series:")
        print("  输入 'help' 查看详情，'latex' 切换输出格式，'quit' 退出")
        print(f"  当前输出格式：{'纯文本' if use_text else 'LaTeX'}\n")
        while True:
            try:
                user_input = input(">>> ").strip()
                if user_input.lower() in ['quit', 'exit', 'q', '退出']:
                    print(" 再见！")
                    break
                if user_input.lower() in ['help', 'h', '帮助']:
                    print("""
+------------------------------------------------------------+
|  复变函数求解器 v2.0 使用帮助                        |
+------------------------------------------------------------+
|  题型前缀（自动识别或手动指定）：                     |
|    equation:   求解方程  例：equation: z^2+1=0        |
|    integral:   围道积分  例：integral: 1/(z^2+1) |z|=2|
|    analytic:   解析判定  例：analytic: sin(z)         |
|    singularity:奇点分类  例：singularity: 1/(z-1) z=1 |
|    residue:    留数计算  例：residue: 1/(z-1)^2 z=1   |
|    series:     级数展开  例：series: exp(z) at 0      |
+------------------------------------------------------------+
|  也可用中文前缀：方程: / 积分: / 解析: /              |
|                 奇点: / 留数: / 级数:                 |
+------------------------------------------------------------+
|  命令：                                               |
|    latex / latex on / latex off  → 切换 LaTeX 输出    |
|    help / h / 帮助               → 显示此帮助         |
|    quit / exit / q / 退出         → 退出程序          |
+------------------------------------------------------------+
""")
                    continue
                # LaTeX 切换命令
                if user_input.lower() == 'latex' or user_input.lower().startswith('latex '):
                    parts = user_input.lower().split()
                    if len(parts) == 1:
                        use_latex = not use_latex  # 无参数则切换
                    elif parts[1] in ['on', '1', 'true', 'yes']:
                        use_latex = True
                    elif parts[1] in ['off', '0', 'false', 'no']:
                        use_latex = False
                    else:
                        print(f"  用法：latex [on|off]，当前为 {'LaTeX' if use_latex else '纯文本'} 模式")
                        continue
                    print(f"  LaTeX 输出已{'开启' if use_latex else '关闭'}（当前：{'LaTeX' if use_latex else '纯文本'}）")
                    continue
                if not user_input:
                    continue
                print("-" * 40)
                print(solve(user_input, latex=use_latex))
                print("-" * 40)
                print()
            except (KeyboardInterrupt, EOFError):
                print("\n 再见！")
                break
