# 🎯 复变函数求解器 — Complex Analysis Solver v2.0

> **考研专用 · 基于 SymPy 符号计算 · 零幻觉设计**

一个基于 Python 符号计算引擎的复变函数求解工具，专为考研/数学系学生设计。所有数学运算由 **SymPy** 执行，**绝不依赖 LLM 编造推导过程**。

## 🆕 v2.0 更新

- **LaTeX 输出**：支持 `--text`（纯文本）、`--markdown`（Markdown）、默认 LaTeX 三种输出模式
- **交互模式**：运行 `python solver.py` 进入交互式 REPL，支持 `help`、`latex on/off` 命令
- **Markdown 格式化**：`--md` / `--markdown` 参数输出带表情符号的 Markdown 结果
- **LaTeX 后处理**：自动将数学表达式包裹 `$$...$$`，清理残余 sympy 格式
- **更完善的题型检测**：新增洛朗级数中文关键词、支点检测、调和性验证
- **围道参数解析增强**：支持 `|z-(a+bi)|=r` 等复数圆心格式
- **方程求解增强**：含共轭/Abs/Re/Im 时自动实虚分离求解

---

## ✨ 功能

| 功能 | 说明 | 示例 |
|:----:|:----|:----|
| 🔢 **方程求解** | 自动因式分解/求根/实虚分离 | `z^6 - 1 = 0` → 6 个单位根 |
| 📐 **围道积分** | 柯西留数定理自动计算 | `∮ 1/(z(1-z)³) dz` |
| ✅ **解析性判定** | C-R 条件逐项验证 | `sin(z)` → 处处解析 ✓ |
| ⚡ **奇点/极点分类** | 自动检测阶数 | `sin(z)/z³` → 2 阶极点 |
| 📍 **留数计算** | 符号留数公式 | `1/(z²+1) at z=i` → -i/2 |
| 📈 **级数展开** | 泰勒/洛朗展开 | `exp(z) at 0 order=6` |

---

## 🚀 快速使用

### 命令行

```bash
# 自动检测题型（默认 LaTeX 输出）
python solver.py "z^6 - 1 = 0"

# 纯文本输出
python solver.py --text "z^2 + 1 = 0"

# Markdown 输出
python solver.py --md "integral: 1/(z(1-z)^3) |z|=0.5"

# 指定题型前缀
python solver.py "integral: 1/(z(1-z)^3) |z|=0.5"
python solver.py "residue: 1/(z^2+1) at z=i"
python solver.py "series: sin(z) at 0 order=8"
python solver.py "singularity: tan(z) at z=pi/2"
python solver.py "analytic: z^2 + 2z"

# 中文前缀同样支持
python solver.py "方程: |z| + z = 1 + i"
python solver.py "积分: 1/(z^2+4) |z|=3"
python solver.py "奇点: 1/(z-1)^2(z+1) at z=1"

# 交互模式（运行后直接输入题目）
python solver.py
```

### 作为 Python 模块

```python
from solver_core import solve

# 方程求解
result = solve("equation: z^6 - 1 = 0")
print(result)

# 围道积分
result = solve("integral: 1/(z^2+1) |z|=2")
print(result)

# 便捷函数
from solver_core import solve_equation, solve_integral, check_analytic
solve_equation("z^2 + 1 = 0")           # 解方程
solve_integral("1/(z^2+1)", 0, 2)       # 积分，圆心0半径2
check_analytic("sin(z)")                 # 判定解析性
```

---

## 📦 安装

```bash
# 依赖
pip install sympy mpmath

# 克隆仓库
git clone https://github.com/lw731946/complex-analysis.git
cd complex-analysis

# 运行
python scripts/solver.py "你的题目"
```

---

## 🧪 示例输出

### 围道积分

```
[STEP] 柯西留数定理计算步骤：
 积分围道：|z - 0.0| = 0.5
 分母零点分析得到极点：[0j, (1+0j)]
 围道内极点及留数：
  z=0j (1阶) → Res = 1+0j
 留数和 = 1
[PASS] ∮f(z)dz = 2πi · ΣRes = 6.283185307*I
```

### 解析性判定

```
[STEP] C-R条件判定步骤：
 f(z) = u(x,y) + i·v(x,y)
 u(x,y) = x**2 - y**2
 v(x,y) = 2*x*y
 ∂u/∂x = 2*x | ∂u/∂y = -2*y
 ∂v/∂x = 2*y | ∂v/∂y = 2*x
 C-R① ∂u/∂x - ∂v/∂y = 0 [PASS]✓
 C-R② ∂u/∂y + ∂v/∂x = 0 [PASS]✓
[PASS] 严格满足C-R条件，函数在复平面上处处解析
```

---

## 🧠 设计哲学

| 原则 | 说明 |
|:----|:-----|
| **零幻觉** | 所有运算由 SymPy 符号计算完成，不会就拒绝，绝不编造答案 |
| **符号优先** | 先尝试解析解、因式分解、留数定理，降级到数值解 |
| **自动验证** | 方程解自动代回原式核验，虚根自动剔除 |
| **安全执行** | 统一异常捕获，不会崩溃 |

---

## 📁 项目结构

```
complex-analysis/
├── SKILL.md               # 技能描述
├── scripts/
│   ├── solver.py          # CLI 入口
│   └── solver_core.py     # 核心求解引擎（SymPy）
├── .gitignore
└── README.md              # 本文档
```

---

## 📄 许可

MIT License
