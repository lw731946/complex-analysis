---
name: "complex-analysis"
description: "复变函数求解器 -- 基于 SymPy 的复变函数计算，支持方程/积分/奇点/留数/级数展开"
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
python solver.py "z^6 - 1 = 0"
python solver.py --text "积分: 1/(z^2+1) |z|=2"
```

## 设计原则

- **零幻觉** — 所有运算由 SymPy 符号计算完成，LLM 不参与数学推理
- **符号优先** — 先尝试解析解、因式分解、留数定理，降级到数值解
- **自动验证** — 方程解自动代回原式核验，虚根自动剔除
- **安全执行** — 统一异常捕获，不会崩溃
- **必须写出完整求解步骤** — 每道题都要给出详细的推导过程

## 调用方式（供 LLM Agent 使用）

```python
import sys, os
sys.path.insert(0, 'scripts/')
from solver_core import solve

result = solve("积分: 1/(z^2+1) |z|=2", latex=False)
print(result)
```

## ⚡ 核心规则：必须写出完整求解步骤

当用户提出复变函数问题时，**必须展示完整的数学推导过程**，不能只给答案。具体做法如下：

### 步骤一：调用 solver 获取结果

将用户题目传入 `solve(题目, latex=False)`，获取得初步结果。

### 步骤二：拆解并写出详细步骤

根据题型，按以下模板写出每一步推导：

#### 方程求解题
```
1. 写出方程：ln(z-1) - cosh(2iz) = ½πi
2. 化简：利用 cosh(2iz) = cos(2z)
3. 整理得：ln(z-1) - cos(2z) = ½πi
4. 两边取指数：z-1 = exp(½πi + cos(2z))
5. 展开复指数：exp(½πi) = i
6. 使用 SymPy 数值求解（nsolve）
```

#### 围道积分题
```
1. 识别被积函数 f(z) 和围道 C
2. 找出围道内的所有奇点
3. 计算各奇点的留数 Res(f, z_k)
4. 由留数定理：∮f(z)dz = 2πi · ΣRes
5. 代入数值计算
```

#### 奇点分类题
```
1. 找出所有使函数不解析的点
2. 对每个奇点 z₀，判断类型：
   - 可去奇点：lim_{z→z₀} f(z) 存在且有限
   - m 阶极点：lim_{z→z₀} (z-z₀)^m·f(z) ≠ 0, ∞
   - 本性奇点：极限不存在且不为 ∞
3. 对极点：求阶数 m 和留数
```

#### 级数展开题
```
1. 确定展开点 z₀ 和展开半径 R
2. 写出泰勒/洛朗级数通项公式
3. 逐项计算系数 a_n
4. 写出前 N 项级数展开式
```

### 步骤三：用 SymPy 辅助计算中间步骤

对于复杂中间计算，直接使用 SymPy 验证：

```python
import sympy as sp

# 验证极限
z = sp.Symbol('z', complex=True)
sp.limit(sp.sin(z)/z, z, 0)  # → 1

# 验证导数
sp.diff(sp.exp(z), z)

# 验证积分
sp.integrate(1/(z**2+1), (z, 0, 1))
```

### 步骤四：给出最终答案

所有步骤完后，用清晰格式给出最终结果。

## 典型完整输出格式

```
【题目】求积分 ∮ 1/(z²+1) dz，其中 |z| = 2

【步骤】
1. 被积函数 f(z) = 1/(z²+1) = 1/[(z-i)(z+i)]
   奇点：z = i, z = -i

2. 围道 |z| = 2 包含两个奇点

3. 计算留数：
   Res(f, i) = lim_{z→i} (z-i)·f(z) = 1/(2i) = -i/2
   Res(f, -i) = lim_{z→-i} (z+i)·f(z) = -1/(2i) = i/2

4. 留数定理：
   ∮ f(z)dz = 2πi · [Res(f, i) + Res(f, -i)]
            = 2πi · (-i/2 + i/2)
            = 2πi · 0
            = 0

【答案】∮ 1/(z²+1) dz = 0

【验证】SymPy 计算确认结果正确 ✅
```

## 文件结构

```
complex-analysis/
  SKILL.md
  scripts/
    solver.py         # CLI 入口
    solver_core.py    # 核心求解引擎
```

## 依赖

- Python 3.8+, sympy, mpmath
- `pip install sympy mpmath`
