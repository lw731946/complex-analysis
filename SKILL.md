---
name: complex-analysis
description: "复变函数求解器 — 解方程、围道积分、解析性判定、奇点/极点分类、留数计算、级数展开(泰勒/洛朗)。基于 SymPy 符号计算，零幻觉设计。"
---

# 复变函数求解器

考研专用复变函数求解工具。基于 SymPy 符号计算引擎，**不会就拒绝，绝不编造答案**。所有代数解自动代回核验。

## 支持题型

| 题型 | 触发词 | 示例 |
|------|--------|------|
| 方程求解 | `equation:` / `方程:` | `equation: z^2 + 1 = 0` |
| 围道积分 | `integral:` / `积分:` | `integral: 1/(z^2+1) \|z\|=2` |
| 解析性判定 | `analytic:` / `解析:` | `analytic: sin(z)` |
| 奇点/极点分类 | `singularity:` / `奇点:` | `singularity: 1/(z-1) z=1` |
| 留数计算 | `residue:` / `留数:` | `residue: 1/(z-1)^2 z=1` |
| 级数展开 | `series:` / `级数:` / `展开:` | `series: exp(z) at 0` |

也支持**无前缀自动识别**，只要有 `=`, `dz`, `弧积分`, `奇点` 等关键词即可。

## 使用方式

```bash
# 直接传入题目（自动检测题型）
python scripts/solver.py "z^6 - 1 = 0"

# 带前缀指定题型
python scripts/solver.py "integral: 1/(z(1-z)^3) |z|=0.5"
python scripts/solver.py "residue: 1/(z^2+1) at z=i"
python scripts/solver.py "series: sin(z) at 0 order=8"

# 中文前缀同样支持
python scripts/solver.py "方程: |z| + z = 1 + i"
python scripts/solver.py "积分: 1/(z^2+4) |z|=3"
python scripts/solver.py "奇点: tan(z) at z=pi/2"

# LaTeX 格式输出（加 --latex 参数）
python scripts/solver.py --latex "z^2 + 1 = 0"
python scripts/solver.py --latex "integral: 1/(z^2+1) |z|=2"
```

## 设计原则

- **零幻觉** — 所有运算由 SymPy 符号计算完成，LLM 不参与数学推理
- **符号优先** — 先尝试解析解、因式分解、留数定理，降级到数值解
- **自动验证** — 方程解自动代回原式核验，虚根自动剔除
- **安全执行** — 统一异常捕获，不会崩溃

## 文件结构

```
complex-analysis/
  SKILL.md
  scripts/
    solver.py         # CLI 入口（供 agent 调用）
    solver_core.py    # 核心求解引擎
```

## 依赖

- Python 3.8+
- sympy
- mpmath

## 典型输出

```
f(z) 的奇点分析
 > z=0 为 1 阶极点
   lim_{z→0} (z)^1·f(z) = 1/2 ≠ 0, ∞

留数计算
[STEP] 一阶极点留数公式：Res = lim(z-z0)f(z)
[PASS] Res(f, 0) = 1/2

围道积分
[STEP] 柯西留数定理计算步骤：
 积分围道：|z| = 2
 围道内极点及留数：
  z=0 (1阶) → Res = 0.5
[PASS] ∮f(z)dz = 2πi · ΣRes = 3.14159265i
```
