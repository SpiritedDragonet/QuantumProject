# coding=utf-8
# cavity.py - 腔系统和光场表示

# 文件权限验证：此注释行确认cavity.py文件已成功设置为可编辑状态

import numpy as np
import qutip as qt
from typing import Dict, List, Optional, Union, Any
import math

class CavitySystem:
    """
    光学腔系统类，用于表示DLCZ协议中的光场和光纤腔
    
    该类实现了腔模式的量子力学表示，支持光子态的创建、变换和测量。
    特别关注DLCZ协议中的信号光子和闲置光子的产生和探测。
    """
    
    def __init__(self, config):
        """
        初始化光学腔系统。
        
        参数：
            config (Config): 配置对象，包含腔参数
        """
        # 从配置中获取腔参数
        self.wavelength = config.get('laser_wavelength', 780e-9)  # 波长 (m)
        self.cavity_length = config.get('cavity_length', 20e-6)   # 腔长 (m)
        self.cavity_waist = config.get('cavity_waist', 2.0e-6)    # 腔模式束腰 (m)
        self.finesse = config.get('cavity_finesse', 1000)         # 腔精细度
        
        # 计算腔参数
        self.free_spectral_range = self._calculate_fsr()          # 自由光谱范围 (Hz)
        self.linewidth = self._calculate_linewidth()              # 腔线宽 (Hz)
        self.kappa = np.pi * self.linewidth                       # 腔场衰减率 (Hz)
        self.mode_volume = self._calculate_mode_volume()          # 模式体积 (m^3)
        
        # 获取原子参数（用于计算耦合强度）
        self.gamma = config.get('gamma_sp', 2.0 * np.pi * 6e6)    # 自发辐射率 (Hz)
        
        # 计算最大耦合强度
        self.g_max = self._calculate_coupling_strength()          # 最大耦合强度 (Hz)
        
        # 光子数量统计
        self.mean_photon_number = 0.0
        self.photon_emission_probability = 0.0
        
        # 创建腔操作符
        self.fock_dim = config.get('fock_dim', 5)                 # 光子数截断维度
        self.a = qt.destroy(self.fock_dim)                        # 湮灭算符
        self.a_dag = self.a.dag()                                 # 产生算符
        self.n = self.a_dag * self.a                              # 数算符
        
        # 设置腔哈密顿量
        cavity_freq = config.get('cavity_frequency', None)
        if cavity_freq is None:
            # 如果未指定腔频率，根据波长计算
            c = 2.99792458e8  # 光速 (m/s)
            self.cavity_frequency = c / self.wavelength  # 腔频率 (Hz)
        else:
            self.cavity_frequency = cavity_freq
        
        self.H_cavity = self.cavity_frequency * self.n            # 腔哈密顿量
        
        # 设置腔跃迁算符
        self.c_ops = [np.sqrt(2 * self.kappa) * self.a]           # 腔损失算符
    
    def _calculate_fsr(self) -> float:
        """
        计算腔的自由光谱范围(FSR)。
        
        FSR = c / (2L)，其中c是光速，L是腔长。
        
        返回值：
            float: 自由光谱范围 (Hz)
        """
        c = 2.99792458e8  # 光速 (m/s)
        return c / (2 * self.cavity_length)
    
    def _calculate_linewidth(self) -> float:
        """
        计算腔线宽。
        
        linewidth = FSR / finesse
        
        返回值：
            float: 腔线宽 (Hz)
        """
        return self.free_spectral_range / self.finesse
    
    def _calculate_mode_volume(self) -> float:
        """
        计算腔模式体积。
        
        V = (pi * w0^2 * L) / 4，其中w0是腔模式束腰半径，L是腔长。
        
        返回值：
            float: 模式体积 (m^3)
        """
        return np.pi * self.cavity_waist**2 * self.cavity_length / 4
    
    def _calculate_coupling_strength(self) -> float:
        """
        计算原子-腔耦合强度。
        
        g = sqrt(3 * lambda^2 * gamma / (8 * pi * V))
        
        其中lambda是波长，gamma是自发辐射率，V是模式体积。
        
        返回值：
            float: 耦合强度 (Hz)
        """
        return np.sqrt(3 * self.wavelength**2 * self.gamma / 
                      (8 * np.pi * self.mode_volume))
    
    def get_position_dependent_coupling(self, position: np.ndarray) -> float:
        """
        计算位置依赖的耦合强度。
        
        g(r) = g_max * exp(-(x^2 + y^2) / w0^2) * cos(k * z)
        
        其中w0是束腰半径，k = 2pi/lambda，z是沿腔轴的位置。
        
        参数：
            position (np.ndarray): 原子位置坐标 [x, y, z] (m)
            
        返回值：
            float: 该位置的耦合强度 (Hz)
        """
        # 提取位置坐标
        x, y, z = position
        
        # 计算横向高斯分布因子
        transverse_factor = np.exp(-(x**2 + y**2) / self.cavity_waist**2)
        
        # 计算轴向驻波因子
        k = 2 * np.pi / self.wavelength
        axial_factor = np.cos(k * z)
        
        # 计算总耦合强度
        g = self.g_max * transverse_factor * axial_factor
        
        return g
    
    def get_atom_cavity_hamiltonian(self, atom_positions: np.ndarray, 
                                   sigma_minus_ops: List[qt.Qobj]) -> qt.Qobj:
        """
        构建原子-腔相互作用哈密顿量。
        
        H_int = sum_i g_i (a^dag * sigma_minus_i + a * sigma_plus_i)
        
        参数：
            atom_positions (np.ndarray): 原子位置数组
            sigma_minus_ops (List[qt.Qobj]): 原子跃迁算符列表
            
        返回值：
            qt.Qobj: 原子-腔相互作用哈密顿量
        """
        H_int = 0
        
        # 为每个原子添加相互作用项
        for i, position in enumerate(atom_positions):
            # 计算该位置的耦合强度
            g_i = self.get_position_dependent_coupling(position)
            
            # 构建相互作用哈密顿量
            # a^dag * sigma_minus + a * sigma_plus
            H_i = g_i * (self.a_dag * sigma_minus_ops[i] + 
                         self.a * sigma_minus_ops[i].dag())
            
            # 累加到总哈密顿量
            H_int += H_i
        
        return H_int
    
    def get_jaynes_cummings_hamiltonian(self, g: float, 
                                      sigma_minus: qt.Qobj) -> qt.Qobj:
        """
        构建Jaynes-Cummings哈密顿量（单原子情况）。
        
        H_JC = g (a^dag * sigma_minus + a * sigma_plus)
        
        参数：
            g (float): 耦合强度
            sigma_minus (qt.Qobj): 原子降算符
            
        返回值：
            qt.Qobj: Jaynes-Cummings相互作用哈密顿量
        """
        return g * (self.a_dag * sigma_minus + self.a * sigma_minus.dag())
    
    def get_tavis_cummings_hamiltonian(self, g_list: List[float], 
                                     sigma_minus_list: List[qt.Qobj]) -> qt.Qobj:
        """
        构建Tavis-Cummings哈密顿量（多原子情况）。
        
        H_TC = sum_i g_i (a^dag * sigma_minus_i + a * sigma_plus_i)
        
        参数：
            g_list (List[float]): 耦合强度列表
            sigma_minus_list (List[qt.Qobj]): 原子降算符列表
            
        返回值：
            qt.Qobj: Tavis-Cummings相互作用哈密顿量
        """
        H_TC = 0
        
        for i, (g_i, sigma_minus_i) in enumerate(zip(g_list, sigma_minus_list)):
            H_i = g_i * (self.a_dag * sigma_minus_i + self.a * sigma_minus_i.dag())
            H_TC += H_i
        
        return H_TC
    
    def get_dlcz_write_hamiltonian(self, g_list: List[float], 
                                  omega_list: List[float],
                                  sigma_eg_list: List[qt.Qobj],
                                  sigma_se_list: List[qt.Qobj]) -> qt.Qobj:
        """
        构建DLCZ协议的写入哈密顿量。
        
        H_write = sum_i [g_i (a^dag * sigma_eg_i) + Omega_i * sigma_se_i] + h.c.
        
        参数：
            g_list (List[float]): 原子-腔耦合强度列表
            omega_list (List[float]): 激光拉比频率列表
            sigma_eg_list (List[qt.Qobj]): |e><g|跃迁算符列表
            sigma_se_list (List[qt.Qobj]): |s><e|跃迁算符列表
            
        返回值：
            qt.Qobj: DLCZ写入哈密顿量
        """
        H_write = 0
        
        for i, (g_i, omega_i, sigma_eg_i, sigma_se_i) in enumerate(
            zip(g_list, omega_list, sigma_eg_list, sigma_se_list)):
            
            # 原子-腔相互作用项
            H_cavity_i = g_i * (self.a_dag * sigma_eg_i + self.a * sigma_eg_i.dag())
            
            # 激光驱动项
            H_laser_i = omega_i * (sigma_se_i + sigma_se_i.dag())
            
            # 累加到总哈密顿量
            H_write += H_cavity_i + H_laser_i
        
        return H_write
    
    def get_dlcz_read_hamiltonian(self, g_list: List[float], 
                                 omega_list: List[float],
                                 sigma_gs_list: List[qt.Qobj],
                                 sigma_es_list: List[qt.Qobj]) -> qt.Qobj:
        """
        构建DLCZ协议的读取哈密顿量。
        
        H_read = sum_i [g_i (a^dag * sigma_gs_i) + Omega_i * sigma_es_i] + h.c.
        
        参数：
            g_list (List[float]): 原子-腔耦合强度列表
            omega_list (List[float]): 激光拉比频率列表
            sigma_gs_list (List[qt.Qobj]): |g><s|跃迁算符列表
            sigma_es_list (List[qt.Qobj]): |e><s|跃迁算符列表
            
        返回值：
            qt.Qobj: DLCZ读取哈密顿量
        """
        H_read = 0
        
        for i, (g_i, omega_i, sigma_gs_i, sigma_es_i) in enumerate(
            zip(g_list, omega_list, sigma_gs_list, sigma_es_list)):
            
            # 原子-腔相互作用项
            H_cavity_i = g_i * (self.a_dag * sigma_gs_i + self.a * sigma_gs_i.dag())
            
            # 激光驱动项
            H_laser_i = omega_i * (sigma_es_i + sigma_es_i.dag())
            
            # 累加到总哈密顿量
            H_read += H_cavity_i + H_laser_i
        
        return H_read
    
    def update_mean_photon_number(self, rho: qt.Qobj) -> float:
        """
        更新并返回腔中的平均光子数。
        
        参数：
            rho (qt.Qobj): 系统密度矩阵或波函数
            
        返回值：
            float: 平均光子数
        """
        # 计算光子数算符的期望值
        self.mean_photon_number = qt.expect(self.n, rho)
        return self.mean_photon_number
    
    def calculate_photon_emission_probability(self, kappa: float, 
                                            dt: float, 
                                            n_photon: float) -> float:
        """
        计算在时间间隔dt内光子从腔中发射的概率。
        
        P_emission = 1 - exp(-2 * kappa * n_photon * dt)
        
        参数：
            kappa (float): 腔场衰减率
            dt (float): 时间间隔
            n_photon (float): 腔内平均光子数
            
        返回值：
            float: 光子发射概率
        """
        self.photon_emission_probability = 1.0 - np.exp(-2 * kappa * n_photon * dt)
        return self.photon_emission_probability
    
    def get_cavity_field_state(self, photon_num: int) -> qt.Qobj:
        """
        获取腔场的量子态。
        
        参数：
            photon_num (int): 光子数量
            
        返回值：
            qt.Qobj: 腔场量子态
        """
        if photon_num == 0:
            return qt.basis(self.fock_dim, 0)  # 真空态
        else:
            return qt.basis(self.fock_dim, photon_num)  # 光子数态
    
    def get_cavity_collapse_operators(self) -> List[qt.Qobj]:
        """
        获取腔的塌缩算符列表。
        
        返回值：
            List[qt.Qobj]: 腔塌缩算符列表
        """
        return self.c_ops
    
    def __str__(self) -> str:
        """返回腔系统的描述字符串"""
        desc = f"光学腔系统\n"
        desc += f"  波长: {self.wavelength*1e9:.2f} nm\n"
        desc += f"  腔长: {self.cavity_length*1e6:.2f} μm\n"
        desc += f"  束腰: {self.cavity_waist*1e6:.2f} μm\n"
        desc += f"  精细度: {self.finesse:.2e}\n"
        desc += f"  自由光谱范围: {self.free_spectral_range/(2*np.pi):.2e} Hz\n"
        desc += f"  线宽: {self.linewidth/(2*np.pi):.2e} Hz\n"
        desc += f"  模式体积: {self.mode_volume:.2e} m^3\n"
        desc += f"  最大耦合强度: {self.g_max/(2*np.pi):.2e} Hz\n"
        desc += f"  平均光子数: {self.mean_photon_number:.2f}\n"
        desc += f"  光子发射概率: {self.photon_emission_probability:.2e}\n"
        return desc 