#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
复变函数求解器 — 考研专用
基于 SymPy 符号计算，零幻觉设计

CLI 调用（供 OpenClaw agent 使用）：
  python solver.py "题目"
  python solver.py "equation: z^2 + 1 = 0"
  python solver.py "integral: 1/(z^2+1) |z|=2"
"""
import sys, os

# 添加当前目录到路径以导入 solver_core
_script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _script_dir)

from solver_core import solve

if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

    if len(sys.argv) > 1:
        problem = ' '.join(sys.argv[1:])
        print(solve(problem))
    else:
        print("用法: python solver.py \"题目\"")
        print("支持: equation: / integral: / analytic: / singularity: / residue: / series:")
