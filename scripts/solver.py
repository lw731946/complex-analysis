#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
复变函数求解器 v2.0 — CLI 入口
==================================
基于 SymPy 符号计算，零幻觉设计

CLI 调用：
  python solver.py "题目"
  python solver.py --text "积分: 1/(z^2+1) |z|=2"
  python solver.py --markdown "singularity: sin(z)/z at z=0"

支持题型前缀：
 equation: / integral: / analytic: / singularity: / residue: / series:
 中文：方程: / 积分: / 解析: / 奇点: / 留数: / 级数: / 展开:
"""

import sys, os

# 添加当前目录到路径以导入 solver_core
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _script_dir)

from solver_core import solve

if __name__ == "__main__":
    # 确保 UTF-8 输出（兼容现代终端）
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

    # 解析 --text / --markdown 标志
    use_text = '--text' in sys.argv
    use_markdown = '--markdown' in sys.argv or '--md' in sys.argv
    args = [a for a in sys.argv[1:] if a not in ['--text', '--markdown', '--md']]
    use_latex = not use_text and not use_markdown

    print("=" * 60)
    print(" 复变函数求解器 v2.0 -- 考研专用")
    mode_label = 'Markdown' if use_markdown else ('纯文本' if use_text else 'LaTeX')
    print(f" 输出模式：{mode_label}")
    print(" 支持：方程 | 积分 | 解析性 | 奇点 | 留数 | 级数展开")
    print("=" * 60)

    if args:
        problem = ' '.join(args)
        print(f"\n[Q] 题目：{problem}")
        print("-" * 40)
        print(solve(problem, latex=use_latex, markdown=use_markdown))
        print("-" * 40)
    else:
        print("\n[STEP] 交互模式")
        print(" 直接输入题目（自动识别题型）")
        print(" 或用前缀指定：equation: / integral: / analytic: /")
        print(" singularity: / residue: / series:")
        print(" 输入 'help' 查看详情，'latex' 切换输出格式，'quit' 退出")
        print(f" 当前输出格式：{'纯文本' if use_text else 'LaTeX'}\n")
        while True:
            try:
                user_input = input(">>> ").strip()
                if user_input.lower() in ['quit', 'exit', 'q', '退出']:
                    print(" 再见！")
                    break
                if user_input.lower() in ['help', 'h', '帮助']:
                    print("""
+------------------------------------------------------------+
| 复变函数求解器 v2.0 使用帮助 |
+------------------------------------------------------------+
| 题型前缀（自动识别或手动指定）： |
| equation: 求解方程 例：equation: z^2+1=0 |
| integral: 围道积分 例：integral: 1/(z^2+1) |z|=2|
| analytic: 解析判定 例：analytic: sin(z) |
| singularity:奇点分类 例：singularity: 1/(z-1) z=1 |
| residue: 留数计算 例：residue: 1/(z-1)^2 z=1 |
| series: 级数展开 例：series: exp(z) at 0 |
+------------------------------------------------------------+
| 也可用中文前缀：方程: / 积分: / 解析: / |
| 奇点: / 留数: / 级数: |
+------------------------------------------------------------+
| 命令： |
| latex / latex on / latex off → 切换 LaTeX 输出 |
| help / h / 帮助 → 显示此帮助 |
| quit / exit / q / 退出 → 退出程序 |
+------------------------------------------------------------+
""")
                    continue
                # LaTeX 切换命令
                if user_input.lower() == 'latex' or user_input.lower().startswith('latex '):
                    parts = user_input.lower().split()
                    if len(parts) == 1:
                        use_latex = not use_latex
                    elif parts[1] in ['on', '1', 'true', 'yes']:
                        use_latex = True
                    elif parts[1] in ['off', '0', 'false', 'no']:
                        use_latex = False
                    else:
                        print(f" 用法：latex [on|off]，当前为 {'LaTeX' if use_latex else '纯文本'} 模式")
                        continue
                    print(f" LaTeX 输出已{'开启' if use_latex else '关闭'}（当前：{'LaTeX' if use_latex else '纯文本'}）")
                    continue
                if not user_input:
                    continue
                print("-" * 40)
                print(solve(user_input, latex=use_latex, markdown=use_markdown))
                print("-" * 40)
                print()
            except (KeyboardInterrupt, EOFError):
                print("\n 再见！")
                break
