# coding=utf-8
# config.py
"""
DLCZ量子接口仿真系统配置模块

本模块实现配置管理，使用单例模式确保全局只有一个配置实例，
允许系统各组件直接访问配置数据而无需传递大量参数。
"""

# 验证文件权限：此注释行确认config.py文件已成功设置为可编辑状态

import os
import numpy as np
from typing import Dict, List, Any, Optional, Union

# 结果存储和调试输出路径
RESULTS_DIR = "results"
DEBUG_OUTPUT_DIR = "docs/debug_output"

class ConfigSingleton:
    """配置单例类，确保全局只有一个配置实例"""
    _instance = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ConfigSingleton, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        if self._initialized:
            return
            
        # 初始化配置
        self._initialize_defaults()
        self._initialized = True
    
    def _initialize_defaults(self):
        """初始化默认配置值"""
        # === 物理常数（所有物理常数均使用SI单位） ===
        self.constants = {
            'hbar': 1.0545718e-34,      # 约化普朗克常数 (J·s)
            'kb': 1.380649e-23,         # 玻尔兹曼常数 (J/K)
            'mu_B': 9.2740100783e-24,   # 玻尔磁子 (J/T)
            'c': 299792458.0,           # 光速 (m/s)
            'epsilon0': 8.8541878128e-12 # 真空介电常数 (F/m)
        }
        
        # === 系统配置 ===
        self.system = {
            'atom_type': 'Rb87',        # 原子类型: 'Rb87' 或 'Cs133'
            'num_atoms': 2,             # 模拟中的原子数量
            'atoms': [],                # 原子参数列表
            'epr_pairs': [],            # EPR对列表
            'simulation_type': 'quantum',  # 模拟类型
            'output_directory': RESULTS_DIR  # 结果输出目录
        }
        
        # === 量子态表示参数 ===
        self.quantum = {
            'levels_per_atom': 3,        # 每个原子的能级数
            'use_full_hilbert_space': True,  # 使用完整希尔伯特空间
            'include_zeeman': False,     # 是否包含Zeeman子态
            'include_rydberg': True,     # 是否包含里德堡态
            'include_lost': True,        # 是否包含丢失态
            'state_representation': 'density_matrix',  # 状态表示方法
            'fock_dim': 3,              # 光场Fock空间截断维度
            'cavity_modes': 1,           # 腔模数量
            'max_photon_number': 5      # 最大光子数
        }
        
        # === 原子参数 (会根据atom_type自动设置) ===
        self.atom = {}  # 将由_load_atom_parameters填充
        
        # === 激光和腔参数 ===
        self.laser = {
            'write_pulse_rabi': 2.0 * np.pi * 5e6,   # 写入脉冲Rabi频率 (Hz)
            'read_pulse_rabi': 2.0 * np.pi * 10e6,   # 读出脉冲Rabi频率 (Hz)
            'write_detuning': 2.0 * np.pi * 10e6,    # 写入失谐 (Hz)
            'read_detuning': 2.0 * np.pi * 10e6,     # 读出失谐 (Hz)
            'write_pulse_duration': 100e-9,          # 写入脉冲持续时间 (s)
            'read_pulse_duration': 100e-9,           # 读出脉冲持续时间 (s)
            'pulse_shape': 'gaussian'                # 脉冲形状
        }
        
        self.cavity = {
            'length': 1e-3,               # 腔长 (m)
            'waist': 30e-6,               # 腔模式束腰 (m)
            'finesse': 1000,              # 腔精细度
            'coupling': 2.0 * np.pi * 1e6,  # 原子-腔耦合强度 (rad/s)
            'kappa': 2.0 * np.pi * 1e5,     # 腔场衰减率 (rad/s)
            'frequency': None              # 腔频率 (Hz)，将自动计算
        }
        
        # === 噪声参数 ===
        self.noise = {
            'gamma_phi': 2.0 * np.pi * 1e3, # 相位弛豫率 (rad/s)
            'gamma_phi_atom': 2.0 * np.pi * 1e3, # 原子相位退相干率 (rad/s)
            'detector_efficiency': 0.5,    # 光子探测器效率
            'detector_dark_count': 10,     # 探测器暗计数率 (Hz)
            'signal_photon_loss': 0.3,     # 信号光子损失率
            'idler_photon_loss': 0.3       # 闲置光子损失率
        }
        
        # === 仿真控制参数 ===
        self.simulation = {
            'total_time': 1e-6,            # 总模拟时间 (s)
            'time_step': 5e-9,             # 时间步长 (s)
            'storage_time': 500e-9,         # 存储时间 (s)
            'use_quantum_jumps': False,     # 是否使用量子跳跃法
            'num_trajectories': 100,        # 量子轨迹数量
            'random_seed': 42              # 随机数种子
        }
        
        # 加载原子参数并初始化
        self._load_atom_parameters()
        self._initialize_atoms()
        self._calculate_cavity_parameters()
        self._ensure_directories_exist()
    
    def _load_atom_parameters(self):
        """加载原子参数"""
        atom_type = self.system['atom_type']
        
        if atom_type == 'Rb87':
            self.atom = {
                'mass': 1.44316060e-25,   # kg
                'dipole_moment': 2.54e-29,  # C·m (D2线)
                'hyperfine_splitting': 6.834682610904e9,  # Hz
                'gamma_sp': 2 * np.pi * 6.065e6,  # Hz (D2线自发辐射率)
                'wavelength_D1': 794.979e-9,  # m
                'wavelength_D2': 780.241e-9,  # m
                'doppler_limit': 146e-6      # K
            }
        elif atom_type == 'Cs133':
            self.atom = {
                'mass': 2.2069517e-25,    # kg
                'dipole_moment': 2.686e-29,  # C·m (D2线)
                'hyperfine_splitting': 9.192631770e9,  # Hz
                'gamma_sp': 2 * np.pi * 5.234e6,  # Hz (D2线自发辐射率)
                'wavelength_D1': 894.593e-9,  # m
                'wavelength_D2': 852.35e-9,   # m
                'doppler_limit': 125e-6       # K
            }
        else:
            raise ValueError(f"不支持的原子类型: {atom_type}")
    
    def _initialize_atoms(self):
        """初始化原子配置"""
        atoms = []
        
        # 添加默认的原子（一对EPR原子）
        # 空间A中的EPR原子
        atom_a = {
            'type': self.system['atom_type'],
            'position': [0, 0, 0],     # 位于空间A原点
            'velocity': np.zeros(3),   # 初始速度为0
            'levels': self.quantum['levels_per_atom'],  # 能级数
            'initial_state': 0,        # 默认在基态 |g⟩
            'region': 'A'              # 所属空间区域
        }
        atoms.append(atom_a)
        
        # 空间B中的EPR原子
        atom_b = {
            'type': self.system['atom_type'],
            'position': [100e-6, 0, 0],  # 位于空间B，距离A 100微米
            'velocity': np.zeros(3),     # 初始速度为0
            'levels': self.quantum['levels_per_atom'],  # 能级数
            'initial_state': 0,          # 默认在基态 |g⟩
            'region': 'B'                # 所属空间区域
        }
        atoms.append(atom_b)
        
        self.system['atoms'] = atoms
        self.system['epr_pairs'] = [(0, 1)]  # 记录默认的EPR对
    
    def _calculate_cavity_parameters(self):
        """计算腔参数"""
        # 只有当腔模式启用时才计算
        if self.quantum['cavity_modes'] > 0:
            # 设置腔频率（如果未设置）
            if self.cavity['frequency'] is None:
                # 默认设置为与原子D2线共振
                self.cavity['frequency'] = self.constants['c'] / self.atom['wavelength_D2']
    
    def _ensure_directories_exist(self):
        """确保必要的目录存在"""
        os.makedirs(self.system['output_directory'], exist_ok=True)
        os.makedirs(DEBUG_OUTPUT_DIR, exist_ok=True)
    
    def get(self, section, key=None, default=None):
        """
        获取配置值，支持点号访问和嵌套访问
        
        参数:
            section (str): 配置区段名称
            key (str, optional): 配置键名，如果为None则返回整个区段
            default: 如果键不存在，返回的默认值
            
        返回:
            配置值或默认值
        """
        if not hasattr(self, section):
            return default
            
        section_dict = getattr(self, section)
        if key is None:
            return section_dict
            
        return section_dict.get(key, default)
    
    def set(self, section, key, value):
        """
        设置配置值
        
        参数:
            section (str): 配置区段名称
            key (str): 配置键名
            value: 配置值
        """
        if not hasattr(self, section):
            setattr(self, section, {})
            
        section_dict = getattr(self, section)
        section_dict[key] = value
        
        # 如果更改了原子类型，重新加载原子参数
        if section == 'system' and key == 'atom_type':
            self._load_atom_parameters()
            
        # 如果更改了腔参数，重新计算派生参数
        if section == 'cavity':
            self._calculate_cavity_parameters()

# 创建全局单例配置实例
config = ConfigSingleton()

def create_default_config():
    """
    获取默认配置实例
    
    返回:
        ConfigSingleton: 全局配置单例
    """
    return config

def get_config_value(section, key=None, default=None):
    """
    获取配置值的便捷函数
    
    参数:
        section (str): 配置区段名称
        key (str, optional): 配置键名，如果为None则返回整个区段
        default: 如果键不存在，返回的默认值
        
    返回:
        配置值或默认值
    """
    return config.get(section, key, default)

def get_physical_constant(name):
    """
    获取物理常数的便捷函数
    
    参数:
        name (str): 常数名称
        
    返回:
        常数值
    """
    return config.constants.get(name)

# 确保目录创建
config._ensure_directories_exist()
