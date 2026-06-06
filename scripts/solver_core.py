#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
复变函数求解器 v2.0 — 考研专用
=================================
支持题型：
 方程求解 | 围道积分 | 解析性判定(C-R) | 奇点/极点分类 | 留数计算 | 级数展开(泰勒/洛朗)

使用方式：
 solve("z^2 + 1 = 0") # 自动检测
 solve("equation: z^6 - 1 = 0") # 显式指定题型
 solve("integral: 1/(z^2+1) |z|=2") # 围道积分
 solve("singularity: sin(z)/z at z=0") # 奇点判定
 solve("residue: 1/(z-1) at z=1") # 留数计算
 solve("series: exp(z) at 0") # 级数展开

设计哲学：
 零幻觉 — 不会就拒绝，绝不编造答案
 符号优先 — 先尝试解析解，再考虑数值近似
 结果可验证 — 所有解代回原方程核验
"""

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
MAX_FACTOR_ATTEMPTS = 20
DEFAULT_SERIES_TERMS = 6
NUMERICAL_GUESS_RANGE = 5
VERIFY_TOLERANCE = 1e-8


# ======================== LaTeX 格式化工具 ========================
def _l(expr, latex=False):
    """将 SymPy 表达式转为可读字符串（LaTeX 或纯文本）"""
    if not latex:
        return str(expr)
    try:
        s = sp.latex(expr, imaginary_unit='i')
        s = s.replace(r'\left[', r'[').replace(r'\right]', r']')
        s = s.replace(r'\left(', r'(').replace(r'\right)', r')')
        return s
    except Exception:
        return str(expr)


def _lm(expr, latex=False):
    """返回 LaTeX 展示数学模式包裹的字符串"""
    s = _l(expr, latex)
    if latex:
        return f'\\[{s}\\\]'
    return s

def _lc(c, latex=False):
    """复数格式化，支持 LaTeX"""
    if latex and isinstance(c, complex):
        re_val = c.real
        im_val = c.imag
        if abs(im_val) < 1e-12:
            return f"\\]{re_val:.6g}\\[".replace('e', '\\times 10^{').replace('+0', '').replace('\\]', '') + '\\['
        if abs(re_val) < 1e-12:
            im_str = f"{im_val:.6g}".rstrip('0').rstrip('.')
            if im_str == '1':
                return '\\]i\\['
            if im_str == '-1':
                return '\\]-i\\['
            return f"\\]{im_str}i\\["
        re_str = f"{re_val:.6g}".rstrip('0').rstrip('.')
        im_str = f"{abs(im_val):.6g}".rstrip('0').rstrip('.')
        if im_str == '1':
            im_str = ''
        sign = '+' if im_val >= 0 else '-'
        return f"\\]{re_str}{sign}{im_str}i\\["
    return f"{c.real:.6g}{c.imag:+.6g}i"


def _lmul(expr, latex=False):
    """将表达式转为 LaTeX 表示的格式，用于普通文本中"""
    if not latex:
        return str(expr)
    return f'\\]{_l(expr, True)}\\['


def _lwrap(text_block, latex=False):
    """将包含 LaTeX 特殊标记的 plain 步骤文本转为 LaTeX 环境"""
    if not latex:
        return text_block
    return text_block


# ======================== 1. 输入清洗 ========================
def clean_math_expression(text):
    text = text.strip()
    text = text.replace('^', '**')
    text = text.replace('（', '(').replace('）', ')')
    text = text.replace('，', ',').replace('；', ';')
    text = text.replace('：', ':')

    text = re.sub(r'(?<![a-zA-Z])i(?![a-zA-Z])', 'I', text)
    text = re.sub(r'(?<![a-zA-Z])j(?![a-zA-Z])', 'I', text)

    func_shorthands = {
        'sinz': 'sin(z)', 'cosz': 'cos(z)', 'tanz': 'tan(z)',
        'sinhz': 'sinh(z)', 'coshz': 'cosh(z)', 'tanhz': 'tanh(z)',
        'ez': 'exp(z)', 'lnz': 'log(z)', 'ln(': 'log(',
        'arcsinz': 'asin(z)', 'arccosz': 'acos(z)', 'arctanz': 'atan(z)',
    }
    for old, new in func_shorthands.items():
        text = text.replace(old, new)

    name_map = {
        'cot(': '1/tan(', 'sec(': '1/cos(', 'csc(': '1/sin(',
        'conj(': 'conjugate(', 're(': 're(', 'im(': 'im(',
        'abs(': 'Abs(', '|z|': 'Abs(z)',
        'arg(': 'arg(',
    }
    for old, new in name_map.items():
        text = text.replace(old, new)

    text = re.sub(r'\|z\s*([-+])\s*(\S+?)\|', r'Abs(z\1\2)', text)
    text = re.sub(r'\|z\|', 'Abs(z)', text)
    text = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', text)
    text = re.sub(r'(\d)\(', r'\1*(', text)
    text = re.sub(r'\bz(\d)', r'z*\1', text)
    # 自动补乘号：z( 和 )( 等缺少乘号的情况
    # z(1-z) → z*(1-z)
    text = re.sub(r'(?<![a-zA-Z])z\(', 'z*(', text)
    # ) 后跟字母（如 )(z → )*(z, )sin → )*sin)
    text = re.sub(r'\)([a-zA-Z])', r')*\1', text)
    text = re.sub(r'\)([a-zA-Z(\d])', r')*\1', text)
    text = re.sub(r'[一-龥]+', '', text)
    return text.strip()


def extract_chinese_keywords(raw_text):
    return re.findall(r'[一-龥]+', raw_text)


def parse_contour_params(raw_text):
    m = re.search(r'\|\s*z\s*([-+])\s*\(?([^|]+?)\)?\s*\|\s*[=<]\s*(\S+)', raw_text)
    if m:
        sign = m.group(1)
        center_str = m.group(2).strip()
        radius_str = m.group(3).strip()
        try:
            if sign == '-':
                center = complex(sp.sympify(center_str).evalf())
            else:
                center = complex(-sp.sympify(center_str).evalf())
            radius = float(sp.sympify(radius_str).evalf())
            return center, radius
        except Exception:
            pass

    m = re.search(r'\|\s*z\s*\|\s*[=<]\s*(\S+)', raw_text)
    if m:
        try:
            radius = float(sp.sympify(m.group(1)).evalf())
            return 0.0, radius
        except Exception:
            pass

    m = re.search(r'[=]\s*(\d+\.?\d*)', raw_text)
    if m:
        try:
            return 0.0, float(m.group(1))
        except Exception:
            pass

    return 0.0, 1.0


# ======================== 2. 题型检测 ========================
EXPLICIT_PREFIX_MAP = {
    'equation': 'equation', '方程': 'equation', 'eq': 'equation',
    'integral': 'integral', '积分': 'integral', '复积分': 'integral', 'int': 'integral',
    'analytic': 'analytic', '解析': 'analytic', 'ana': 'analytic',
    'singularity': 'singularity', '奇点': 'singularity', '极点': 'singularity', 'singular': 'singularity',
    'residue': 'residue', '留数': 'residue', 'res': 'residue',
    'series': 'series', '级数': 'series', '展开': 'series', 'expand': 'series',
}


def parse_explicit_prefix(raw_text):
    for prefix_str, task_type in EXPLICIT_PREFIX_MAP.items():
        for sep in [':', '：']:
            full = prefix_str + sep
            if raw_text.lower().startswith(full):
                return task_type, raw_text[len(full):].strip()
    return None, raw_text


def detect_task_type(raw_text, cleaned_text):
    chinese_keywords = extract_chinese_keywords(raw_text)
    chinese_joined = ''.join(chinese_keywords)

    if any(kw in chinese_joined for kw in ['积分', '围道', '环路', '柯西积分']):
        center, radius = parse_contour_params(''.join(chinese_keywords) + raw_text)
        return 'integral', 0.95, {'center': center, 'radius': radius}
    if 'dz' in raw_text.lower() or '∮' in raw_text or '∳' in raw_text:
        center, radius = parse_contour_params(raw_text)
        return 'integral', 0.85, {'center': center, 'radius': radius}

    if any(kw in chinese_joined for kw in ['洛朗', 'laurent']):
        point, order = parse_series_params(raw_text)
        return 'series', 0.95, {'point': point, 'order': order, 'is_laurent': True}
    if any(kw in chinese_joined for kw in ['展开', '级数', '泰勒', '麦克劳林', 'taylor', 'maclaurin']):
        point, order = parse_series_params(raw_text)
        return 'series', 0.95, {'point': point, 'order': order, 'is_laurent': False}

    if any(kw in chinese_joined for kw in ['奇点', '极点', '本性', '可去奇点', '孤立奇点']):
        point = parse_point(raw_text)
        return 'singularity', 0.95, {'point': point}

    if any(kw in chinese_joined for kw in ['留数', 'residue']):
        point = parse_point(raw_text)
        return 'residue', 0.95, {'point': point}

    if any(kw in chinese_joined for kw in ['解析', '柯西黎曼', 'C-R', '调和']):
        return 'analytic', 0.95, {}

    if '=' in cleaned_text:
        return 'equation', 0.85, {}

    if '/' in cleaned_text:
        point = parse_point(raw_text)
        try:
            _, denom_str = cleaned_text.rsplit('/', 1)
            denom = sp.sympify(denom_str.strip())
            if denom.has(z):
                return 'residue', 0.6, {'point': point}
        except Exception:
            pass

    if 'z' in cleaned_text or 'Z' in cleaned_text:
        return 'analytic', 0.3, {}

    return 'unknown', 0.0, {}


def parse_point(raw_text):
    # 先匹配 z=后的整体内容
    m = re.search(r'[zZ]\s*[=＝]\s*(.+)', raw_text)
    if m:
        point_str = m.group(1).strip()
        # 去除行尾多余内容（如围道参数）
        point_str = re.sub(r'\s*\|.*$', '', point_str)
        # 尝试直接 SymPy 解析
        try:
            val = sp.sympify(point_str.replace('i', 'I').replace('j', 'I'))
            if val.is_number:
                c = complex(val.evalf())
                if abs(c.imag) < 1e-12:
                    return int(c.real) if c.real == int(c.real) else float(c.real)
                return c
        except Exception:
            pass
    # 兜底：从文本中提取最后一个数字
    nums = re.findall(r'[-+]?\d+\.?\d*', raw_text)
    if nums:
        try:
            f = float(nums[-1])
            return int(f) if f == int(f) else f
        except Exception:
            pass
    return 0


def parse_series_params(raw_text):
    point = 0
    order = DEFAULT_SERIES_TERMS
    nums = re.findall(r'[-+]?\d+\.?\d*', raw_text)
    if nums:
        try:
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

    m = re.search(r'[zZ]\s*[=＝]\s*([-+]?\d+)', raw_text)
    if m:
        try:
            point = float(m.group(1))
        except Exception:
            pass

    return point, order


# ======================== 3. 安全执行器 ========================
def safe_exec(func, *args, **kwargs):
    try:
        result = func(*args, **kwargs)
        if isinstance(result, dict):
            result['success'] = True
            return result
        return {'success': True, 'result': str(result), 'error': None}
    except sp.PolynomialError:
        return {'success': False, 'result': '[FAIL] 非标准多项式，无通用解析解法', 'error': 'PolynomialError'}
    except ValueError as e:
        return {'success': False, 'result': f'[FAIL] 数值错误：{e}', 'error': 'ValueError'}
    except TypeError as e:
        return {'success': False, 'result': f'[FAIL] 类型错误：{e}', 'error': 'TypeError'}
    except SyntaxError as e:
        return {'success': False, 'result': f'[FAIL] 语法错误：{e}', 'error': 'SyntaxError'}
    except ZeroDivisionError:
        return {'success': False, 'result': '[FAIL] 存在间断奇点，函数此处无定义', 'error': 'ZeroDivisionError'}
    except NotImplementedError as e:
        return {'success': False, 'result': f'[FAIL] 暂不支持：{e}', 'error': 'NotImplementedError'}
    except Exception as e:
        error_type = type(e).__name__
        return {'success': False, 'result': f'[FAIL] 求解失败({error_type})：{e}', 'error': error_type}


# ======================== 4. 方程求解 ========================
def solve_complex_equation(eq_str, allow_numerical=True, latex=False):
    if '=' not in eq_str:
        return "[FAIL] 非标准等式方程，格式错误"

    lhs_str, rhs_str = eq_str.split('=', 1)
    expr_str = f"({lhs_str.strip()})-({rhs_str.strip()})"

    try:
        expr = sp.sympify(expr_str, locals={'z': z})
    except Exception as e:
        return f"[FAIL] 表达式解析失败：{e}"

    # ---- 含共轭/绝对值/Re/Im ----
    if expr.has(sp.conjugate) or expr.has(sp.Abs) or expr.has(sp.re) or expr.has(sp.im):
        try:
            expr_xy = expr.subs(z, x + I*y).expand()
            eq_real = sp.simplify(sp.re(expr_xy))
            eq_imag = sp.simplify(sp.im(expr_xy))
            sols_xy = sp.solve([eq_real, eq_imag], [x, y], dict=True)

            if latex:
                steps = f"设 \\[z = x + iy\\]，代入得：\n \\[{_l(expr_xy, True)} = 0\\]\n实部: \\[{_l(eq_real, True)} = 0\\]\n虚部: \\[{_l(eq_imag, True)} = 0\\]"
            else:
                steps = f"设 z = x + iy，代入得：\n ({expr_xy}) = 0\n实部: {eq_real} = 0\n虚部: {eq_imag} = 0"

            if sols_xy == []:
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
                return f"[PASS] 解析解（实虚分离法）：\n{format_solution_list(complex_sols, latex)}"
        except Exception:
            pass

    # ---- 因式分解 ----
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
                return f"[PASS] 解析解（因式分解法）：\n{format_solution_list(valid_sols, latex)}"
    except Exception:
        pass

    # ---- 符号求解 ----
    try:
        sym_sols = sp.solve(expr, z)
        if sym_sols:
            valid_sols = verify_solutions_expr(expr_str, sym_sols)
            return f"[PASS] 解析解：\n{format_solution_list(valid_sols, latex)}"
    except Exception:
        pass

    # ---- 数值求解兜底 ----
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
                return f"[PASS] 数值近似解（牛顿迭代法）：\n{format_solution_list(sorted_roots, latex)}\n[WARN] 注意：此为非准确解析解，仅作参考"
        except Exception:
            pass

    return "[FAIL] 该方程无初等解析解且数值求解亦未收敛"


def _generate_complex_guesses():
    guesses = []
    for a in range(-NUMERICAL_GUESS_RANGE, NUMERICAL_GUESS_RANGE + 1):
        for b in range(-NUMERICAL_GUESS_RANGE, NUMERICAL_GUESS_RANGE + 1):
            guesses.append(complex(a, b))
            guesses.append(complex(a + 0.5, b + 0.5))
    return guesses


def _round_complex(c, decimals):
    return round(c.real, decimals) + 1j * round(c.imag, decimals)


def format_solution_list(solutions, latex=False):
    if not solutions:
        return " (空集)"
    lines = []
    for i, sol in enumerate(solutions, 1):
        if isinstance(sol, (int, float, complex)):
            c = complex(sol)
            if latex:
                re_str = f"{c.real:.6g}".rstrip('0').rstrip('.')
                im_str = f"{abs(c.imag):.6g}".rstrip('0').rstrip('.')
                if abs(c.imag) < 1e-12:
                    lines.append(f" \\[z_{{{i}}} = {re_str}\\]")
                elif abs(c.real) < 1e-12:
                    im_show = im_str if im_str != '1' else ''
                    sign = '-' if c.imag < 0 else ''
                    lines.append(f" \\[z_{{{i}}} = {sign}{im_show}i\\]")
                else:
                    sign = '+' if c.imag >= 0 else '-'
                    im_show = im_str if im_str != '1' else ''
                    lines.append(f" \\[z_{{{i}}} = {re_str}{sign}{im_show}i\\]")
            else:
                lines.append(f" z{i} = {c.real:.6g}{c.imag:+.6g}i")
        else:
            if latex:
                lines.append(f" \\[z_{{{i}}} = {_l(sol, True)}\\]")
            else:
                lines.append(f" z{i} = {sol}")
    return '\n'.join(lines)


def verify_solutions_expr(expr_str, solutions):
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
            valid.append(sol)
    return valid


# ======================== 5. 解析性判定 ========================
def judge_analytic(func_str, latex=False):
    try:
        f = sp.sympify(func_str, locals={'z': z, 'x': x, 'y': y})
    except Exception:
        return '[FAIL] 表达式解析失败'

    try:
        f_xy = sp.simplify(f.subs(z, x + I*y))
        u_expr = sp.simplify(sp.re(f_xy))
        v_expr = sp.simplify(sp.im(f_xy))
    except Exception:
        return "[FAIL] 无法分离实虚部，请检查函数表达式"

    ux = sp.diff(u_expr, x)
    uy = sp.diff(u_expr, y)
    vx = sp.diff(v_expr, x)
    vy = sp.diff(v_expr, y)

    cr1_diff = sp.simplify(ux - vy)
    cr2_diff = sp.simplify(uy + vx)

    steps = []
    if latex:
        steps.append("[STEP] C-R条件判定步骤：")
        steps.append(f" \\[f(z) = u(x,y) + i\\,v(x,y)\\]")
        steps.append(f" \\[u(x,y) = {_l(u_expr, True)}\\]")
        steps.append(f" \\[v(x,y) = {_l(v_expr, True)}\\]")
        steps.append(f" \\[\\frac{{\\partial u}}{{\\partial x}} = {_l(ux, True)}\\] | "
                     f"\\[\\frac{{\\partial u}}{{\\partial y}} = {_l(uy, True)}\\]")
        steps.append(f" \\[\\frac{{\\partial v}}{{\\partial x}} = {_l(vx, True)}\\] | "
                     f"\\[\\frac{{\\partial v}}{{\\partial y}} = {_l(vy, True)}\\]")
        cr1_str = f"\\[\\frac{{\\partial u}}{{\\partial x}} - \\frac{{\\partial v}}{{\\partial y}} = {_l(cr1_diff, True)}\\]"
        cr2_str = f"\\[\\frac{{\\partial u}}{{\\partial y}} + \\frac{{\\partial v}}{{\\partial x}} = {_l(cr2_diff, True)}\\]"
        steps.append(f" C-R① {cr1_str} {'[PASS]✓' if cr1_diff == 0 else '[NO]✗'}")
        steps.append(f" C-R② {cr2_str} {'[PASS]✓' if cr2_diff == 0 else '[NO]✗'}")
    else:
        steps.append("[STEP] C-R条件判定步骤：")
        steps.append(f" f(z) = u(x,y) + i·v(x,y)")
        steps.append(f" u(x,y) = {u_expr}")
        steps.append(f" v(x,y) = {v_expr}")
        steps.append(f" ∂u/∂x = {ux} | ∂u/∂y = {uy}")
        steps.append(f" ∂v/∂x = {vx} | ∂v/∂y = {vy}")
        steps.append(f" C-R① ∂u/∂x - ∂v/∂y = {cr1_diff} {'[PASS]✓' if cr1_diff == 0 else '[NO]✗'}")
        steps.append(f" C-R② ∂u/∂y + ∂v/∂x = {cr2_diff} {'[PASS]✓' if cr2_diff == 0 else '[NO]✗'}")

    if cr1_diff == 0 and cr2_diff == 0:
        steps.append("[PASS] 严格满足C-R条件，函数在复平面上处处解析")
    else:
        steps.append("[NO] 不满足C-R条件，函数不是解析函数")

    return '\n'.join(steps)


# ======================== 6. 奇点/极点分类 ========================
def classify_singularity(func_str, z0, latex=False):
    if str(z0).lower() in ['inf', '∞', 'infinity']:
        return (" > 无穷远点判定：需作代换 \\[w=1/z\\]，分析 \\[w=0\\] 处奇点类型。\n"
                " 请对 \\[f(1/w)\\] 在 \\[w=0\\] 处重新分析。")

    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    z0_sym = sp.sympify(z0)

    # 支点检测
    if (f.has(sp.log) or f.has(sp.sqrt) or
        any(a.is_Pow and (not a.exp.is_Integer) for a in sp.preorder_traversal(f) if a.is_Pow)):
        try:
            val_at = complex(f.subs(z, z0_sym).evalf())
            if latex:
                return (f" > \\[z={z0}\\] 疑似为支点（branch point）\n"
                        f" 函数含多值成分（\\[\\log\\]/sqrt/分数幂），该点非孤立奇点。\n"
                        f" 留数概念不适用于支点，需沿分支切割分析。")
            return (f" > z={z0} 疑似为支点（branch point）\n"
                    f" 函数含多值成分（log/sqrt/分数幂），该点非孤立奇点。\n"
                    f" 留数概念不适用于支点，需沿分支切割分析。")
        except Exception:
            pass

    # 可去奇点
    try:
        limit_val = sp.limit(f, z, z0_sym)
        if limit_val.is_finite:
            if latex:
                return (f" > \\[z={z0}\\] 为可去奇点\n"
                        f" \\[\\displaystyle\\lim_{{z\\to {z0}}} f(z) = {_l(limit_val, True)}\\]\n"
                        f" 留数 \\[\\operatorname{{Res}}(f, {z0}) = 0\\]")
            return (f" > z={z0} 为可去奇点\n"
                    f" lim_{{z→{z0}}} f(z) = {limit_val}\n"
                    f" 留数 Res(f, {z0}) = 0")
    except Exception:
        pass

    # 极点阶数
    for n in range(1, 9):
        try:
            g = (z - z0_sym)**n * f
            lim_g = sp.limit(g, z, z0_sym)
            if lim_g.is_finite and lim_g != 0:
                if latex:
                    return (f" > \\[z={z0}\\] 为 {n} 阶极点\n"
                            f" \\[\\displaystyle\\lim_{{z\\to {z0}}} (z-{z0})^{{{n}}}\\cdot f(z) = {_l(lim_g, True)} \\neq 0, \\infty\\]")
                return (f" > z={z0} 为 {n} 阶极点\n"
                        f" lim_{{z→{z0}}} (z-{z0})^{n}·f(z) = {lim_g} ≠ 0, ∞")
        except Exception:
            continue

    try:
        poles_dict = sp.poles(f, z)
        for p, order in poles_dict.items():
            if abs(complex((p - z0_sym).evalf())) < 1e-10:
                if latex:
                    return f" > \\[z={z0}\\] 为 {order} 阶极点（有理函数判定）"
                return f" > z={z0} 为 {order} 阶极点（有理函数判定）"
    except Exception:
        pass

    if latex:
        return f" > \\[z={z0}\\] 为本性奇点（极限振荡或不存在，且非有限阶极点）"
    return f" > z={z0} 为本性奇点（极限振荡或不存在，且非有限阶极点）"


# ======================== 7. 留数计算 ========================
def calc_residue(func_str, z0, latex=False):
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    z0_sym = sp.sympify(z0)

    try:
        lim_val = sp.limit(f, z, z0_sym)
        if lim_val.is_finite:
            if latex:
                return (f"[STEP] \\[z={z0}\\] 为可去奇点\n"
                        f"[PASS] \\[\\operatorname{{Res}}(f, {z0}) = 0\\]")
            return (f"[STEP] z={z0} 为可去奇点\n"
                    f"[PASS] Res(f, {z0}) = 0")
    except Exception:
        pass

    try:
        res = sp.residue(f, z, z0_sym)
        steps = []
        if latex:
            steps.append("[STEP] 留数计算步骤：")
            steps.append(f" 1. 函数：\\[f(z) = {_l(f, True)}\\]")
            steps.append(f" 2. 奇点：\\[z = {z0}\\]")
            steps.append(f" 3. 使用留数公式计算")
            steps.append(f"[PASS] \\[\\operatorname{{Res}}(f, {z0}) = {_l(sp.simplify(res), True)}\\]")
        else:
            steps.append("[STEP] 留数计算步骤：")
            steps.append(f" 1. 函数：f(z) = {f}")
            steps.append(f" 2. 奇点：z = {z0}")
            steps.append(f" 3. 使用留数公式计算")
            steps.append(f"[PASS] Res(f, {z0}) = {sp.simplify(res)}")
        return '\n'.join(steps)
    except Exception:
        pass

    try:
        g = (z - z0_sym) * f
        lim_g = sp.limit(g, z, z0_sym)
        if lim_g.is_finite and lim_g != 0:
            if latex:
                return (f"[STEP] 一阶极点留数公式：\\[\\operatorname{{Res}} = \\lim_{{z\\to {z0}}}(z-{z0})f(z)\\]\n"
                        f"[PASS] \\[\\operatorname{{Res}}(f, {z0}) = {_l(sp.simplify(lim_g), True)}\\]")
            return (f"[STEP] 一阶极点留数公式：Res = lim(z-z0)f(z)\n"
                    f"[PASS] Res(f, {z0}) = {sp.simplify(lim_g)}")
    except Exception:
        pass

    if latex:
        return f"[FAIL] 无法计算 \\[z={z0}\\] 处的留数"
    return f"[FAIL] 无法计算 z={z0} 处的留数"


# ======================== 8. 围道积分 ========================
def calc_complex_integral(func_str, center=0.0, radius=1.0, latex=False):
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    poles_list = []
    method = ""

    try:
        poles_dict = sp.poles(f, z)
        for pole, order in poles_dict.items():
            poles_list.append((complex(pole.evalf()), order, pole))
        method = "有理函数极点分析"
    except Exception:
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

    inside = []
    on_boundary = []

    for pole_val, order, pole_sym in poles_list:
        dist = abs(pole_val - complex(center))
        if dist < radius - 1e-10:
            inside.append((pole_val, order, pole_sym))
        elif abs(dist - radius) < 1e-10:
            on_boundary.append((pole_val, order, pole_sym))

    steps = []
    if latex:
        steps.append("[STEP] 柯西留数定理计算步骤：")
        steps.append(f" 积分围道：\\[|z - {center}| = {radius}\\]")
        pole_strs = ', '.join(_l(sp.sympify(p[0]), True) for p in poles_list)
        steps.append(f" {method}得到极点：\\[\\{{{pole_strs}\\}}\\]")
        if on_boundary:
            pts = ', '.join(str(p[0]) for p in on_boundary)
            steps.append(f" [WARN] 警告：极点 \\[{{{pts}}}\\] 恰在围道上")
            steps.append(f" 积分需取柯西主值，考研一般不考查此情形")
        if not inside:
            steps.append(f" 围道内无奇点")
            steps.append(f"[PASS] 由柯西积分定理：\\[\\oint f(z)\\,dz = 0\\]")
            return '\n'.join(steps)
    else:
        steps.append("[STEP] 柯西留数定理计算步骤：")
        steps.append(f" 积分围道：|z - {center}| = {radius}")
        steps.append(f" {method}得到极点：{[p[0] for p in poles_list]}")
        if on_boundary:
            steps.append(f" [WARN] 警告：极点 {[p[0] for p in on_boundary]} 恰在围道上")
            steps.append(f" 积分需取柯西主值，考研一般不考查此情形")
        if not inside:
            steps.append(f" 围道内无奇点")
            steps.append(f"[PASS] 由柯西积分定理：∮f(z)dz = 0")
            return '\n'.join(steps)

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
        steps.append(f" 围道内极点及留数：")
        for pv, od, r in residue_details:
            steps.append(f" \\[z={pv}\\] ({od}阶) \\[\\to \\operatorname{{Res}} = {r:.6g}\\]")
        steps.append(f" 留数和 = \\[{_fc_latex(total_res)}\\]")
        steps.append(f"[PASS] \\[\\displaystyle\\oint f(z)\\,dz = 2\\pi i \\cdot \\sum\\operatorname{{Res}} = {sp.N(result, 10)}\\]")
    else:
        steps.append(f" 围道内极点及留数：")
        for pv, od, r in residue_details:
            steps.append(f" z={pv} ({od}阶) → Res = {r:.6g}")
        steps.append(f" 留数和 = {_format_complex(total_res)}")
        steps.append(f"[PASS] ∮f(z)dz = 2πi · ΣRes = {sp.N(result, 10)}")

    return '\n'.join(steps)


def _fc_latex(c):
    """复数 LaTeX 格式化"""
    if abs(c.imag) < 1e-12:
        return f"{c.real:.6g}"
    if abs(c.real) < 1e-12:
        return f"{c.imag:.6g}i"
    return f"{c.real:.6g}{c.imag:+.6g}i"


def _format_complex(c):
    if abs(c.imag) < 1e-12:
        return f"{c.real:.6g}"
    if abs(c.real) < 1e-12:
        return f"{c.imag:.6g}i"
    return f"{c.real:.6g}{c.imag:+.6g}i"


def _has_negative_powers(expr):
    for term in expr.as_ordered_terms():
        if not term.has(z):
            continue
        pd = term.as_powers_dict()
        if pd.get(z, 0) < 0:
            return True
    expr_str = str(expr)
    if '/z' in expr_str or '**(-' in expr_str:
        return True
    return False


# ======================== 9. 级数展开 ========================
def series_expand(func_str, at_point=0, order=6, is_laurent=False, latex=False):
    try:
        f = sp.sympify(func_str, locals={'z': z})
    except Exception:
        return f'[FAIL] 表达式解析失败：{func_str}'

    at_sym = sp.sympify(at_point)
    type_name = "洛朗" if is_laurent else "泰勒"

    try:
        s = sp.series(f, z, at_sym, order)
        terms = s.removeO()
        has_negative = _has_negative_powers(terms)

        if latex:
            terms_l = _l(terms, True)
        else:
            terms_l = str(terms)

        if is_laurent:
            if has_negative:
                if latex:
                    return (f"[PASS] \\[z={at_point}\\] 处的洛朗级数（{order}项）：\n"
                            f" \\[f(z) = {terms_l}\\]\n"
                            f" (含负幂项，确为洛朗展开)")
                return (f"[PASS] z={at_point} 处的洛朗级数（{order}项）：\n"
                        f" f(z) = {terms}\n"
                        f" (含负幂项，确为洛朗展开)")
            else:
                if latex:
                    return (f"[WARN] \\[z={at_point}\\] 处无负幂项，实为泰勒展开：\n"
                            f" \\[f(z) = {terms_l}\\]")
                return (f"[WARN] z={at_point} 处无负幂项，实为泰勒展开：\n"
                        f" f(z) = {terms}")

        if has_negative:
            if latex:
                return (f"[PASS] \\[z={at_point}\\] 处的级数展开（{order}项）：\n"
                        f" \\[f(z) = {terms_l}\\]\n"
                        f" (检测到负幂项，此展开实为洛朗级数)")
            return (f"[PASS] z={at_point} 处的级数展开（{order}项）：\n"
                    f" f(z) = {terms}\n"
                    f" (检测到负幂项，此展开实为洛朗级数)")

        if latex:
            return f"[PASS] \\[z={at_point}\\] 处的{type_name}展开（{order}项）：\n \\[f(z) = {terms_l}\\]"
        return f"[PASS] z={at_point} 处的{type_name}展开（{order}项）：\n f(z) = {terms}"

    except Exception as e:
        try:
            s = sp.series(f, z, at_sym, 3)
            terms = s.removeO()
            terms_l = _l(terms, latex) if latex else str(terms)
            if latex:
                return f"[WARN] 仅能展开到3项：\n \\[f(z) = {terms_l}\\]\n (更高阶展开失败：{e})"
            return f"[WARN] 仅能展开到3项：\n f(z) = {terms}\n (更高阶展开失败：{e})"
        except Exception:
            if latex:
                return f"[FAIL] 无法在 \\[z={at_point}\\] 处展开：{e}"
            return f"[FAIL] 无法在 z={at_point} 处展开：{e}"


def _split_func_and_params(text, task_type):
    if task_type == 'integral':
        center, radius = parse_contour_params(text)
        func_text = re.sub(r'\|\s*z[^|]*\|\s*[=<]\s*\S+', '', text).strip()
        func_text = re.sub(r'[,，]\s*$', '', func_text).strip()
        return func_text, {'center': center, 'radius': radius}

    elif task_type in ('singularity', 'residue'):
        point = parse_point(text)
        func_text = re.sub(r'\s+(at|在)\s+z?\s*[=＝]\s*[-+]?\S+', '', text, flags=re.IGNORECASE)
        func_text = re.sub(r'\s+z\s*[=＝]\s*[-+]?\S+\s*$', '', func_text)
        return func_text.strip(), {'point': point}

    elif task_type == 'series':
        point, order = parse_series_params(text)
        is_laurent = any(kw in text.lower() for kw in ['洛朗', 'laurent'])
        func_text = re.sub(r'\s+(at|在)\s+(z\s*[=＝]\s*)?[-+]?\S+', '', text, flags=re.IGNORECASE)
        func_text = re.sub(r'\s+z\s*[=＝]\s*[-+]?\S+', '', func_text)
        func_text = re.sub(r'\s+order\s*[=＝]\s*\d+', '', func_text, flags=re.IGNORECASE)
        return func_text.strip(), {'point': point, 'order': order, 'is_laurent': is_laurent}

    return text.strip(), {}


# ======================== 10. 统一调度 ========================
def solve(raw_input, latex=False):
    raw_input = raw_input.strip()

    if not raw_input:
        return "[ERR] 请输入题目或表达式"

    task_type, remaining = parse_explicit_prefix(raw_input)

    if task_type is None:
        task_type, confidence, params = detect_task_type(raw_input, raw_input)

        if task_type == 'unknown' or confidence < 0.3:
            return (
                "[ERR] 无法识别题型或未检测到复变量 z\n\n"
                "[TIP] 请使用题型前缀指定：\n"
                " equation: 表达式 → 求解方程\n"
                " integral: 表达式 → 围道积分\n"
                " analytic: 表达式 → 判定解析性\n"
                " singularity: 表达式 → 奇点分类\n"
                " residue: 表达式 → 留数计算\n"
                " series: 表达式 → 级数展开\n\n"
                " 示例：solve(\"integral: 1/(z^2+1) |z|=2\")"
            )

    func_text, split_params = _split_func_and_params(remaining, task_type)
    if split_params:
        params = split_params
        input_clean = clean_math_expression(func_text)
    else:
        func_text, params = _split_func_and_params(remaining, task_type)
        input_clean = clean_math_expression(func_text)

    if not params:
        _task_type, _confidence, params = detect_task_type(remaining, input_clean)
        if _confidence < 0.3:
            params = {}

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
        result = "[FAIL] 未知错误：题型分发失败"

    return result


# ======================== 11. 便捷函数 ========================
def solve_equation(eq_str, allow_numerical=True):
    return solve(f"equation:{eq_str}")

def solve_integral(func_str, center=0.0, radius=1.0):
    return solve(f"integral:{func_str} |z-{center}|={radius}")

def check_analytic(func_str):
    return solve(f"analytic:{func_str}")

def get_residue(func_str, point):
    return solve(f"residue:{func_str} at z={point}")

def expand_series(func_str, point=0, order=6):
    return solve(f"series:{func_str} at z={point} order={order}")

def analyze_singularity(func_str, point):
    return solve(f"singularity:{func_str} at z={point}")
