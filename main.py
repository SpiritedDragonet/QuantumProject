"""
中性原子量子接口仿真系统主模块

本模块作为量子接口仿真系统的入口点，负责创建和管理仿真进程与UI进程。
它仅提供基本的命令行参数解析和进程管理功能，具体的仿真逻辑和UI实现
分别在utils/run_simulation.py和utils/ui.py中实现。

文件权限验证：此注释行用于验证文件已成功设置为可编辑状态。
"""

# coding=utf-8
import argparse
import multiprocessing as mp
import os
import time
import numpy as np
import logging

# 导入配置模块
from config import create_default_config

# 导入UI模块
from utils.ui import run_visualization

# 导入核心模块，用于调试模拟
from core.quantum_state import QuantumState, Frame
from core.time_evolution import TimeEvolution
from core.system_state import initialize_quantum_systems
from core.lindblad_solver import LindbladSolver

# 导入run_debug_simulation_with_frames（从utils.run_simulation导入以避免重复代码）
from utils.run_simulation import run_debug_simulation_with_frames

def main():
    """
    主程序入口点，负责解析命令行参数并创建相应进程
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='中性原子量子接口仿真系统')
    parser.add_argument('--nogui', action='store_true', help='仅执行仿真，不启动UI界面')
    parser.add_argument('--load', type=str, help='加载指定的仿真结果进行可视化', default=None)
    args = parser.parse_args()
    
    # 获取配置
    config = create_default_config()
    
    # 场景1: 仅加载结果进行可视化
    if args.load:
        print(f"从{args.load}加载仿真结果并启动UI")
        run_visualization(config, load_path=args.load)
        return
    
    # 导入run_simulation（延迟导入避免循环导入问题）
    from utils.run_simulation import run_simulation
    
    # 场景2: 仅执行仿真（无UI）
    if args.nogui:
        print("开始执行仿真（无可视化）")
        # 不传递debug_frames参数，避免循环导入
        run_simulation(config)
        return
    
    # 场景3: 同时执行仿真和可视化（默认）
    print("同时执行仿真和可视化")
    
    # 创建进程间通信队列
    data_queue = mp.Queue()
    
    # 创建并启动仿真和可视化进程
    sim_process = mp.Process(target=run_simulation, args=(config, data_queue))
    vis_process = mp.Process(target=run_visualization, args=(config, data_queue))
    
    sim_process.start()
    vis_process.start()
    
    # 等待仿真进程完成
    sim_process.join()
    print("仿真进程已完成")
    
    # 提示用户关闭UI窗口
    if vis_process.is_alive():
        print("可视化进程仍在运行，请手动关闭UI窗口")

if __name__ == "__main__":
    main()
