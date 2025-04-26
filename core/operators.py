# coding=utf-8
# noise_models.py

"""
噪声模型与量子退相干模块 (基于 QuTiP)

本模块负责根据物理参数构建 Lindblad 崩塌算符 (c_ops) 列表，
用于描述量子系统的各种退相干和耗散过程。

支持新的 QuantumState 结构，可以为多种不同大小的量子系统构建噪声模型。

函数调用关系图：
```mermaid
graph TD
    %% LindbladSolver类
    solver_init[LindbladSolver.__init__] --> None
    solver_evolve[LindbladSolver.evolve_density_matrix] --> None
    solver_solve[LindbladSolver.solve] --> None
    solver_stats[LindbladSolver.get_performance_stats] --> None
    
    %% NoiseModelBuilder类
    init[NoiseModelBuilder.__init__] --> None
    build_collapse[NoiseModelBuilder.build_collapse_operators] --> add_spont_em
    build_collapse --> add_cavity
    
    %% 算符添加方法
    add_spont_em[NoiseModelBuilder.add_spontaneous_emission] --> _add_op
    add_cavity[NoiseModelBuilder.add_cavity_decay] --> _add_op
    add_dephasing[NoiseModelBuilder.add_dephasing] --> _add_op
    _add_op[NoiseModelBuilder._add_op_to_list] --> get_operator[QuantumState.get_operator]
    
    %% 静态工厂方法
    factory_build[NoiseModelBuilder.build_for_system] --> build_collapse
    factory_multi[NoiseModelBuilder.build_multi_system_noise] --> factory_build
```
"""

import numpy as np
import qutip as qt
import time
import logging
from typing import Dict, List, Tuple, Any, Optional, Union, Callable, Set

# 导入配置模块和量子态模块
from config import get_config_value
from core.quantum_state import QuantumState

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LindbladSolver:
    """
    Lindblad主方程求解器
    
    实现量子开放系统的Lindblad主方程求解，包括相干演化和退相干过程。
    该类主要使用QuTiP的求解器，并添加了特定于DLCZ协议的功能。
    """
    
    def __init__(self):
        """
        初始化Lindblad求解器
        
        被调用者: 外部调用, TimeEvolution.__init__
        调用: 无其他函数
        """
        # 性能统计
        self.solve_time = 0
        self.solve_count = 0
        
    def evolve_density_matrix(self, rho: qt.Qobj, hamiltonian: qt.Qobj, 
                           c_ops: List[qt.Qobj], dt: float) -> qt.Qobj:
        """
        单步演化密度矩阵
        
        使用QuTiP的演化函数进行单步Lindblad演化
        
        被调用者: TimeEvolution.evolve_system, TimeEvolution.evolve_all_systems
        调用: qt.mesolve
        
        参数:
            rho: 当前密度矩阵
            hamiltonian: 哈密顿量
            c_ops: 崩塌算符列表
            dt: 时间步长
            
        返回:
            qt.Qobj: 演化后的密度矩阵
        """
        start_time = time.time()
        
        # 只计算一步，从t=0到t=dt
        tlist = np.array([0, dt])
        
        # 设置选项：禁用进度条，不存储所有状态，只需最终状态
        options = {
            'store_states': False,
            'progress_bar': False
        }
        
        # 使用QuTiP的mesolve求解
        result = qt.mesolve(hamiltonian, rho, tlist, c_ops, [], options=options)
        
        # 更新统计信息
        self.solve_time += time.time() - start_time
        self.solve_count += 1
        
        # 返回最终状态
        return result.states[-1]
        
    def solve(self, hamiltonian: qt.Qobj, rho0: qt.Qobj, tlist: np.ndarray, 
             c_ops: List[qt.Qobj], e_ops: Optional[List[qt.Qobj]] = None) -> Any:
        """
        求解Lindblad主方程
        
        使用QuTiP的mesolve求解器，计算给定哈密顿量和跃迁算符下系统的演化。
        
        被调用者: TimeEvolution.evolve_multiple_steps
        调用: qt.mesolve
        
        参数:
            hamiltonian: 系统哈密顿量
            rho0: 初始量子态（密度矩阵或态矢量）
            tlist: 时间点数组
            c_ops: 跃迁算符列表
            e_ops: 期望值算符列表
            
        返回:
            qt.solver.Result: 求解结果，包含随时间演化的态和期望值
        """
        start_time = time.time()
        
        # 设置默认选项
        options = {
            'store_states': True,
            'atol': 1e-8,
            'rtol': 1e-6
        }
        
        # 使用QuTiP的mesolve求解
        result = qt.mesolve(hamiltonian, rho0, tlist, c_ops, e_ops, options=options)
        
        # 更新统计信息
        self.solve_time += time.time() - start_time
        self.solve_count += 1
        
        return result
        
    def get_performance_stats(self) -> Dict[str, float]:
        """
        获取求解器性能统计信息
        
        被调用者: 外部调用
        调用: 无其他函数
        
        返回:
            Dict[str, float]: 包含性能统计的字典
        """
        avg_time = self.solve_time / max(1, self.solve_count)
        
        return {
            'total_solve_time': self.solve_time,
            'solve_count': self.solve_count,
            'avg_solve_time': avg_time
        }

class NoiseModelBuilder:
    """
    构建 Lindblad 崩塌算符 (c_ops) 列表。
    
    根据配置参数（如衰减率、退相干率）和量子态信息，
    组装描述各种噪声过程的 QuTiP 算符列表。
    支持多种系统尺寸和类型。
    """
    
    def __init__(self, quantum_state: QuantumState):
        """
        初始化噪声模型构建器。
        
        被调用者: 外部调用, build_for_system
        调用: 无其他函数
        
        参数:
            quantum_state: 关联的 QuantumState 对象，提供维度和算符构建支持。
            
        抛出:
            TypeError: 如果quantum_state不是QuantumState实例
        """
        if not isinstance(quantum_state, QuantumState):
            raise TypeError("必须提供一个 QuantumState 实例")
        self.qs = quantum_state
        self.c_ops: List[qt.Qobj] = []
        
        logger.info(f"NoiseModelBuilder 初始化完成，关联到系统 '{self.qs.system_name}'")

    def _add_op_to_list(self, base_op: qt.Qobj, rate: float, target_subsystem_idx: int):
        """
        将局域算符乘以 sqrt(rate)，扩展到全局空间，并添加到 c_ops 列表。
        
        被调用者: add_spontaneous_emission, add_cavity_decay, add_dephasing
        调用: QuantumState.get_operator
        
        参数:
            base_op: 局部算符
            rate: 噪声速率
            target_subsystem_idx: 目标子系统的索引
        """
        if rate > 1e-15: # 只添加速率大于阈值的算符
            global_op = self.qs.get_operator(np.sqrt(rate) * base_op, target_subsystem_idx)
            self.c_ops.append(global_op)
            logger.debug(f"添加崩塌算符: target={self.qs.subsystem_names[target_subsystem_idx]}, rate={rate:.4g}")

    def add_spontaneous_emission(self, decay_rate: float, atom_idx: int, level_trans: Tuple[int, int]):
        """
        添加原子的自发辐射通道。
        
        被调用者: build_collapse_operators
        调用: _add_op_to_list
        
        参数:
            decay_rate: 衰减率 (gamma)
            atom_idx: 目标原子的索引
            level_trans: 跃迁的能级对 (level_final, level_initial)，例如 (0, 2) 代表 |0⟩⟨2|
        """
        level_final, level_initial = level_trans
        # 检查能级是否有效
        atom_dim = self.qs.dims[atom_idx]
        if not (0 <= level_final < atom_dim and 0 <= level_initial < atom_dim):
             logger.warning(f"原子 {atom_idx} 的自发辐射跃迁能级 {level_trans} 无效 (原子维度 {atom_dim})，跳过此通道。")
             return
             
        # 构建局域跃迁算符 |final⟩⟨initial|
        local_op = qt.basis(atom_dim, level_final) * qt.basis(atom_dim, level_initial).dag()
        # 添加到 c_ops 列表 (自动乘以 sqrt(rate) 并扩展)
        self._add_op_to_list(local_op, decay_rate, atom_idx)
        logger.info(f"添加原子 {atom_idx} 自发辐射: |{level_final}⟩⟨{level_initial}|, rate={decay_rate:.4g}")

    def add_cavity_decay(self, kappa: float, cavity_idx: int):
        """
        添加腔场的衰减通道。
        
        被调用者: build_collapse_operators
        调用: _add_op_to_list
        
        参数:
            kappa: 腔场衰减率
            cavity_idx: 目标腔场的索引
        """
        # 检查腔维度
        cavity_dim = self.qs.dims[cavity_idx]
        # 构建腔场湮灭算符
        a = qt.destroy(cavity_dim)
        # 添加到 c_ops 列表
        self._add_op_to_list(a, kappa, cavity_idx)
        logger.info(f"添加腔 {cavity_idx} 衰减: 速率 kappa={kappa:.4g}")

    def add_dephasing(self, dephasing_rate: float, target_subsystem_idx: int, op_type: str = 'sz'):
        """
        添加退相干通道，包括相位退相干和去极化噪声。
        
        被调用者: build_collapse_operators
        调用: _add_op_to_list
        
        参数:
            dephasing_rate: 退相干率
            target_subsystem_idx: 目标子系统的索引
            op_type: 算符类型 ('sz'=相位退相干, 'sx/sy/sz'=去极化)
        """
        target_dim = self.qs.dims[target_subsystem_idx]
        
        if op_type == 'sz':
            # 对角算符，例如原子的 σz 或腔场的光子数算符
            if target_dim >= 2:
                if self.qs.subsystem_types[target_subsystem_idx] == 'atom':
                    # 对于原子，使用 σz 算符
                    diag_elements = np.zeros(target_dim)
                    for i in range(target_dim):
                        diag_elements[i] = 1 if i % 2 == 1 else -1
                else:
                    # 对于其他类型，使用编号算符
                    diag_elements = np.arange(target_dim)
                    
                op = qt.Qobj(np.diag(diag_elements))
                self._add_op_to_list(op, dephasing_rate, target_subsystem_idx)
                logger.info(f"添加子系统 {target_subsystem_idx} 的相位退相干: rate={dephasing_rate:.4g}")
                
        elif op_type in ['sx', 'sy'] and target_dim >= 2:
            # σx 或 σy 算符
            if op_type == 'sx':
                op = qt.sigmax() if target_dim == 2 else qt.jmat((target_dim-1)/2, 'x')
            else:
                op = qt.sigmay() if target_dim == 2 else qt.jmat((target_dim-1)/2, 'y')
                
            self._add_op_to_list(op, dephasing_rate, target_subsystem_idx)
            logger.info(f"添加子系统 {target_subsystem_idx} 的 σ{op_type[1]} 退相干: rate={dephasing_rate:.4g}")

    def build_collapse_operators(self) -> List[qt.Qobj]:
        """
        构建所有噪声通道的崩塌算符列表。
        
        被调用者: 外部调用, build_for_system
        调用: add_spontaneous_emission, add_cavity_decay, add_dephasing,
              QuantumState.get_subsystems_by_type
        
        返回:
            List[qt.Qobj]: 所有崩塌算符列表
        """
        # 清空现有的算符列表
        self.c_ops = []
        
        # 添加原子噪声
        atom_indices = self.qs.get_subsystems_by_type('atom')
        if atom_indices:
            # 从配置中获取原子噪声参数
            gamma_sp = get_config_value('atom', 'gamma_sp', default=2.0 * np.pi * 6e6)  # 自发辐射率 (rad/s)
            gamma_dephasing = get_config_value('atom', 'gamma_dephasing', default=2.0 * np.pi * 1e3)  # 退相干率 (rad/s)
            
            # 添加原子自发辐射和退相干
            for atom_idx in atom_indices:
                atom_dim = self.qs.dims[atom_idx]
                
                # 添加主要的自发辐射通道
                if atom_dim >= 3: # 至少有 |g₁⟩, |g₂⟩, |e₁⟩
                    self.add_spontaneous_emission(gamma_sp, atom_idx, (0, 2))     # |e₁⟩ → |g₁⟩
                    self.add_spontaneous_emission(gamma_sp * 0.5, atom_idx, (1, 2))  # |e₁⟩ → |g₂⟩
                
                # 如果有第二个激发态，添加相应的通道
                if atom_dim >= 4: # 有 |e₂⟩
                    self.add_spontaneous_emission(gamma_sp * 0.5, atom_idx, (0, 3))  # |e₂⟩ → |g₁⟩
                    self.add_spontaneous_emission(gamma_sp, atom_idx, (1, 3))     # |e₂⟩ → |g₂⟩
                    
                # 如果有里德堡态，添加其衰减通道
                if atom_dim >= 5: # 有 |r⟩
                    for i in range(2):  # 衰减到两个基态
                        self.add_spontaneous_emission(gamma_sp * 0.1, atom_idx, (i, 4))
                    for i in range(2, 4):  # 衰减到两个激发态
                        self.add_spontaneous_emission(gamma_sp * 0.05, atom_idx, (i, 4))
                        
                # 为所有能级添加退相干
                self.add_dephasing(gamma_dephasing, atom_idx, 'sz')
        
        # 添加腔噪声
        cavity_indices = self.qs.get_subsystems_by_type('cavity')
        if cavity_indices:
            # 从配置中获取腔噪声参数
            kappa = get_config_value('cavity', 'kappa', default=2.0 * np.pi * 1e6)  # 腔衰减率 (rad/s)
            
            # 添加腔光子损失算符
            for cavity_idx in cavity_indices:
                self.add_cavity_decay(kappa, cavity_idx)
        
        logger.info(f"为系统 '{self.qs.system_name}' 构建了 {len(self.c_ops)} 个崩塌算符")
        return self.c_ops
    
    @staticmethod
    def build_for_system(quantum_state: QuantumState, 
                        noise_types: Optional[List[str]] = None) -> List[qt.Qobj]:
        """
        为给定的量子系统构建噪声模型。
        
        被调用者: 外部调用, TimeEvolution.evolve_with_system_hamiltonian,
                TimeEvolution.evolve_multiple_steps_with_hamiltonian,
                build_multi_system_noise
        调用: NoiseModelBuilder.__init__, NoiseModelBuilder.build_collapse_operators
        
        参数:
            quantum_state: 量子系统
            noise_types: 要包含的噪声类型列表，如为None则包含所有类型
            
        返回:
            List[qt.Qobj]: 崩塌算符列表
        """
        # 使用构建器创建基本噪声列表
        builder = NoiseModelBuilder(quantum_state)
        all_cops = builder.build_collapse_operators()
        
        # 如果未指定噪声类型或要求包含所有类型，直接返回完整列表
        if noise_types is None:
            return all_cops
            
        # 否则，只保留特定类型的噪声
        filtered_cops = []
        
        # 根据系统中的子系统和指定的噪声类型进行过滤
        have_atoms = any(t == 'atom' for t in quantum_state.subsystem_types)
        have_cavities = any(t == 'cavity' for t in quantum_state.subsystem_types)
        
        # 对每个崩塌算符应用过滤器
        for cop in all_cops:
            if ('spontaneous_emission' in noise_types and have_atoms) or \
               ('dephasing' in noise_types and have_atoms) or \
               ('cavity_loss' in noise_types and have_cavities) or \
               ('depolarizing' in noise_types):
                filtered_cops.append(cop)
        
        # 如果过滤后为空，但应该有噪声，返回原始列表
        if not filtered_cops and all_cops:
            logger.warning(f"噪声类型过滤后为空，返回所有噪声: {noise_types}")
            return all_cops
            
        return filtered_cops
    
    @staticmethod
    def build_multi_system_noise(systems: Dict[str, QuantumState], 
                                noise_types: Optional[List[str]] = None) -> Dict[str, List[qt.Qobj]]:
        """
        为多个量子系统构建噪声模型。
        
        被调用者: 外部调用, TimeEvolution.evolve_all_systems_with_hamiltonian
        调用: build_for_system
        
        参数:
            systems: 量子系统字典 {系统名称: 量子系统}
            noise_types: 要包含的噪声类型列表，如果为None则包含所有类型
            
        返回:
            Dict[str, List[qt.Qobj]]: 噪声算符字典 {系统名称: 崩塌算符列表}
        """
        noise_dict = {}
        for name, system in systems.items():
            noise_dict[name] = NoiseModelBuilder.build_for_system(system, noise_types)
        return noise_dict 