# coding=utf-8
# run_simulation.py

"""
仿真运行器模块

本模块负责执行量子接口仿真的主要逻辑，包括初始化量子系统、
执行DLCZ协议和量子隐形传态、收集结果等。它是仿真进程的核心组件。
"""

# 导入必要的库
import time
import sys
import os
import traceback
from multiprocessing import Queue
import logging
import numpy as np
from typing import Dict, Any, List, Optional

# 导入配置模块
from config import config

# 导入核心组件
from core.quantum_state import QuantumState, Frame
from core.time_evolution import TimeEvolution
from core.system_state import SystemState, initialize_quantum_systems
from core.lindblad_solver import LindbladSolver
from core.entanglement import calculate_concurrence

# 导入物理过程模块
from physics.atomic_process import AtomicProcess
from physics.photon_process import PhotonProcess

# 导入协议模块
from protocols.initialization import initialize_protocol
from protocols.two_cavity_epr import implement_dlcz_protocol
from protocols.quantum_teleport import run_teleportation_protocol as run_teleportation

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("仿真运行器")

def run_simulation(config: Dict[str, Any], data_queue=None) -> str:
    """
    执行量子接口仿真的主要入口函数
    
    参数:
        config: 配置字典
        data_queue: 可选的进程间通信队列，用于向可视化进程发送数据
    
    返回:
        str: 仿真结果的保存路径
    """
    simulation_start_time = time.time()
    logger.info("开始执行量子接口仿真...")
    
    # 确保结果目录存在
    results_dir = "results"
    os.makedirs(results_dir, exist_ok=True)
    
    try:
        # 初始化量子系统状态
        logger.info("初始化量子系统...")
        physical_params = {
            'atom_dim': 5,  # 5能级原子
            'cavity_dim': 3,  # 3能级腔场
            'cavity_coupling': 25.0,  # 耦合强度 (MHz)
            'cavity_kappa': 6.0,  # 腔衰减率 (MHz)
            'atomic_gamma': 6.066,  # 原子自发辐射率 (MHz)
            'dephasing_rate': 0.1  # 退相干率 (MHz)
        }
        system_B, system_C = initialize_quantum_systems(physical_params)
        
        # 创建时间演化控制器
        controller_config = {
            'save_directory': results_dir,
            'simulation_id': f'sim_{int(time.time())}',
            'max_frames_in_memory': 500,
            'frames_per_save': 100
        }
        solver = LindbladSolver()
        evolution = TimeEvolution(solver=solver, config=controller_config)
        
        # 添加系统到演化控制器
        evolution.add_system(system_B)
        evolution.add_system(system_C)
        logger.info(f"已添加两个系统到演化控制器: {system_B.system_name}, {system_C.system_name}")
        
        # 设置DLCZ参数
        dlcz_params = {
            'rabi_freq_B': 2.0 * np.pi,  # 拉比频率 (MHz)
            'rabi_freq_C': 2.0 * np.pi,
            'detuning_B': 0.0,  # 失谐 (MHz)
            'detuning_C': 0.0,
            'pulse_shape': 'gaussian',
            'pulse_duration': 50.0,  # 脉冲持续时间 (ns)
            'pulse_delay': 10.0,  # 脉冲延迟 (ns)
            'pulse_width': 15.0,  # 脉冲宽度 (ns)
            'phase_B': 0.0,  # 相位 (rad)
            'phase_C': 0.0
        }
        
        # 执行DLCZ协议
        logger.info("开始执行DLCZ协议产生远程纠缠...")
        dlcz_result = evolution.run_epr_generation_process(
            simulation_time=200.0e-9,  # 200纳秒
            num_steps=100,
            pulse_params=dlcz_params,
            collect_frames=True,
            max_attempts=3,  # 最多尝试3次
            system_name_B="B",
            system_name_C="C"
        )
        
        # 处理DLCZ结果
        if dlcz_result['success']:
            logger.info(f"DLCZ协议成功: 波包重叠度 γ = {dlcz_result.get('gamma', 'N/A')}")
            
            # 如果有数据队列，向可视化进程发送数据
            if data_queue is not None:
                data_queue.put({
                    'event_type': 'dlcz_completed',
                    'success': True,
                    'gamma': dlcz_result.get('gamma', 0),
                    'frames_path': evolution.get_current_frames_path()
                })
            
            # 执行量子隐形传态
            logger.info("开始执行量子隐形传态...")
            teleport_params = {
                'initial_state_A': 'superposition',  # 或者 'ground', 'excited'
                'bell_measurement_type': 'ideal'
            }
            
            teleport_result = run_teleportation(
                state_manager=evolution,
                params=teleport_params
            )
            
            if teleport_result['success']:
                logger.info(f"量子隐形传态成功: 保真度 = {teleport_result.get('fidelity', 'N/A')}")
                
                # 如果有数据队列，向可视化进程发送数据
                if data_queue is not None:
                    data_queue.put({
                        'event_type': 'teleport_completed',
                        'success': True,
                        'fidelity': teleport_result.get('fidelity', 0),
                        'frames_path': evolution.get_current_frames_path()
                    })
            else:
                logger.warning("量子隐形传态失败")
                if data_queue is not None:
                    data_queue.put({
                        'event_type': 'teleport_completed',
                        'success': False,
                        'error': teleport_result.get('error', 'Unknown error')
                    })
        else:
            logger.warning(f"DLCZ协议失败: {dlcz_result.get('error', 'Unknown error')}")
            if data_queue is not None:
                data_queue.put({
                    'event_type': 'dlcz_completed',
                    'success': False,
                    'error': dlcz_result.get('error', 'Unknown error')
                })
        
        # 保存最终结果
        final_result_path = os.path.join(results_dir, f"{controller_config['simulation_id']}_final_result.npz")
        np.savez(
            final_result_path,
            dlcz_success=dlcz_result.get('success', False),
            gamma=dlcz_result.get('gamma', 0),
            teleport_success=teleport_result.get('success', False) if 'teleport_result' in locals() else False,
            fidelity=teleport_result.get('fidelity', 0) if 'teleport_result' in locals() else 0,
            simulation_time=time.time() - simulation_start_time
        )
        
        logger.info(f"仿真完成，总用时: {time.time() - simulation_start_time:.2f}秒")
        logger.info(f"结果已保存到: {final_result_path}")
        
        # 返回结果目录
        return results_dir
        
    except Exception as e:
        logger.error(f"仿真执行出错: {e}")
        logger.error(traceback.format_exc())
        
        # 通知UI进程发生错误
        if data_queue is not None:
            data_queue.put({
                'event_type': 'simulation_error',
                'error': str(e),
                'traceback': traceback.format_exc()
            })
        
        return results_dir

def run_debug_simulation_with_frames(config, frames_count=None, duration_ns=None, is_self_check=False):
    """
    运行生成帧数据的调试模拟
    
    参数:
        config: 配置字典
        frames_count: 要生成的帧数，如果为None使用默认值2000
        duration_ns: 模拟总时长（纳秒），如果为None使用默认值1000
        is_self_check: 是否作为自检的一部分运行，如果是则简化部分过程
    
    返回:
        str: 保存的结果文件路径
    """
    prefix = "自检" if is_self_check else "调试模拟"
    frames_count = frames_count or 2000
    duration_ns = duration_ns or 1000.0
    
    print("="*60)
    print(f"开始{prefix}，生成{frames_count}帧时间演化数据（时长：{duration_ns}ns）...")
    print("="*60)
    
    # 确保结果目录存在
    results_dir = "results/debug"
    if is_self_check:
        results_dir = "results/self_check"
    os.makedirs(results_dir, exist_ok=True)
    
    # 创建时间演化控制器 - 优化配置以高效保存更多帧
    controller_config = {
        'save_directory': results_dir,
        'simulation_id': f'{"self_check" if is_self_check else "debug"}_simulation_{int(time.time())}',
        'max_frames_in_memory': min(1000, frames_count),  # 减少内存占用
        'frames_per_save': min(1000, frames_count)        # 调整保存频率
    }
    solver = LindbladSolver()
    evolution = TimeEvolution(solver=solver, config=controller_config)
    
    # 初始化量子系统
    print(f"{prefix}：初始化量子系统...")
    physical_params = {
        'atom_dim': 5,  # 5能级原子
        'cavity_dim': 3,  # 3能级腔场
        'cavity_coupling': 25.0,  # 耦合强度 (MHz)
        'cavity_kappa': 6.0,  # 腔衰减率 (MHz)
        'atomic_gamma': 6.066,  # 原子自发辐射率 (MHz)
        'dephasing_rate': 0.1  # 退相干率 (MHz)
    }
    system_B, system_C = initialize_quantum_systems(physical_params)
    
    # 添加系统到演化控制器
    evolution.add_system(system_B)
    evolution.add_system(system_C)
    print(f"添加了两个系统到演化控制器: {system_B.system_name}, {system_C.system_name}")
    
    # 设置模拟参数
    total_time = duration_ns * 1e-9  # 转换为秒
    num_frames = frames_count
    dt = total_time / num_frames
    
    print(f"开始执行时间演化，总时间: {total_time*1e9:.1f} ns, 帧数: {num_frames}...")
    print(f"帧间隔: {dt*1e9:.3f} ns")
    
    # 模拟参数 - 适合Rb原子激发的参数
    pulse_params = {
        'rabi_freq_B': 2.0 * np.pi,  # 拉比频率 (MHz)
        'rabi_freq_C': 2.0 * np.pi,
        'detuning_B': 0.0,  # 失谐 (MHz)
        'detuning_C': 0.0,
        'pulse_shape': 'gaussian',
        'pulse_duration': 50.0,  # 脉冲持续时间 (ns)
        'pulse_delay': 10.0,  # 脉冲延迟 (ns)
        'pulse_width': 15.0,  # 脉冲宽度 (ns)
        'phase_B': 0.0,  # 相位 (rad)
        'phase_C': 0.0
    }
    
    try:
        # 执行EPR纠缠生成过程
        result = evolution.run_epr_generation_process(
            system_name_B="B",
            system_name_C="C",
            simulation_time=total_time,
            num_steps=num_frames,
            pulse_params=pulse_params,
            collect_frames=True,
            max_attempts=1  # 仅尝试一次，专注于生成帧
        )
        
        # 确保保存所有剩余帧
        if len(evolution.frames) > 0:
            print(f"保存剩余的 {len(evolution.frames)} 帧...")
            evolution._save_frames()
        
        # 获取保存的文件路径
        frames_path = evolution.get_current_frames_path()
        
        # 保存结果数据
        result_file = os.path.join(results_dir, f"{controller_config['simulation_id']}_result.json")
        import json
        with open(result_file, 'w') as f:
            # 将不可序列化对象转换为字符串
            serializable_result = {}
            for key, value in result.items():
                if key == 'final_state' or key == 'frames':
                    serializable_result[key] = str(value)
                elif isinstance(value, np.ndarray):
                    serializable_result[key] = value.tolist()
                else:
                    serializable_result[key] = value
            
            # 添加帧保存路径信息
            serializable_result['frames_path'] = frames_path
            
            json.dump(serializable_result, f, indent=2)
        
        # 打印结果统计
        frame_count = evolution.frame_count
        print("="*60)
        print(f"{prefix}完成！生成了 {frame_count} 帧")
        print(f"结果文件已保存到: {result_file}")
        if 'wave_functions' in result:
            print(f"波包重叠度 γ = {result.get('gamma', 'N/A')}")
        print(f"帧文件保存在: {frames_path}")
        print("="*60)
        
        # 如果不是自检，添加额外的验证步骤
        if not is_self_check:
            print("验证保存的帧数据...")
            from utils.visualization import verify_frame_files
            
            verification_result = verify_frame_files(
                directory_path=results_dir,
                load_frames=True,
                plot_sample=False
            )
            
            if verification_result['success']:
                print(f"验证成功! 找到 {verification_result['valid_files']} 个有效帧文件，包含 {verification_result['frames_count']} 帧数据")
            else:
                print(f"验证失败: {verification_result.get('error', '未知错误')}")
        
        return results_dir
    
    except Exception as e:
        import traceback
        print(f"执行时间演化时出错: {e}")
        traceback.print_exc()
        print(f"尽管出错，仍尝试保存当前帧")
        
        # 保存已收集的帧
        if hasattr(evolution, 'frames') and evolution.frames:
            evolution._save_frames()
            print(f"已保存 {len(evolution.frames)} 帧到 {results_dir}")
        
        return results_dir

def main():
    """
    主函数，当作为独立脚本运行时调用
    """
    # 运行仿真
    results_dir = run_simulation(config)
    print(f"仿真结果已保存在: {results_dir}")
    return results_dir

if __name__ == "__main__":
    main()

# 导出当前配置，供其他模块使用
config = config