"""
量子态表示与希尔伯特空间管理模块

本模块负责维护量子系统的状态表示和希尔伯特空间的管理，包括：
- 希尔伯特空间的初始化与维护
- 密度矩阵的设置与获取
- 物理参数的注册与获取

函数调用关系：
```mermaid
graph TD
    A[QuantumState.__init__] --> B[init_hilbert_space]
    C[set_density_matrix] --> D[检查密度矩阵合法性]
    F[register_parameter] --> G[添加物理参数]
    H[get_parameter] --> I[获取物理参数]
    N[get_density_matrix] --> O[获取密度矩阵]
```

数学背景：
- 全局希尔伯特空间 $\mathcal{H} = \mathcal{H}_1 \otimes \mathcal{H}_2 \otimes \cdots \otimes \mathcal{H}_n$
- 密度矩阵 $\rho \in \mathbb{C}^{D \times D}$，其中 $D = \prod_i d_i$

系统描述：
- 原子系统：5个态（两个初态、两个激发态、一个里德堡态）
- 腔系统：3个态（|0⟩, |1⟩, |2⟩ 光子数态）
- 复合系统组合：
  * 单原子系统：5×5密度矩阵
  * 原子-腔系统：5×3=15维（15×15密度矩阵）
  * 双原子系统：5×5=25维（25×25密度矩阵）
  * 三原子系统：5×5×5=125维（125×125密度矩阵）
"""

import numpy as np
import scipy.linalg as la
from typing import List, Dict, Tuple, Optional, Union, Callable, Any
import warnings
import enum
import qutip as qt
import time
import datetime

# 定义系统类型枚举（可选使用）
class SubsystemType(enum.Enum):
    """子系统类型枚举，用于标识子系统的物理含义"""
    ATOM = "atom"        # 原子系统
    CAVITY = "cavity"    # 腔系统
    OTHER = "other"      # 其他系统类型

class QuantumState:
    """
    量子态表示类

    该类负责管理量子系统的状态表示，包括希尔伯特空间的结构、密度矩阵、
    以及相关的物理参数。

    主要属性：
        subsystem_dims (List[int]): 子系统的维度列表
        subsystem_types (List[SubsystemType]): 子系统的类型列表（可选）
        total_dim (int): 总的希尔伯特空间维度
        rho (qt.Qobj): 密度矩阵，QuTiP的Qobj对象
        hamiltonian (qt.Qobj): 哈密顿量矩阵，QuTiP的Qobj对象
        lindblad_ops (List[qt.Qobj]): Lindblad算符列表，每个是QuTiP的Qobj对象
        parameters (Dict): 存储物理参数的字典
        state_name (str): 量子态的名称描述
    """

    def __init__(self, state_name: str = "default"):
        """
        初始化QuantumState对象

        初始化一个空的量子态对象，需要通过init_hilbert_space方法
        设置希尔伯特空间结构后才能使用。

        参数：
            state_name (str): 量子态名称，默认为"default"
        """
        self.subsystem_dims = None  # 子系统维度列表，如[5, 3]表示两个子系统，维度分别为5和3
        self.subsystem_types = None  # 子系统类型列表，如[ATOM, CAVITY]，可选
        self.total_dim = None       # 总希尔伯特空间维度 D_total = prod(subsystem_dims)
        self.rho = None             # 密度矩阵 ρ，QuTiP的Qobj对象
        self.hamiltonian = None     # 哈密顿量 H，QuTiP的Qobj对象
        self.lindblad_ops = None    # Lindblad算符列表，每个是QuTiP的Qobj对象
        self.parameters = {}        # 存储物理参数的字典，如激光参数、衰减率等
        self.subsystem_cumulative_dims = None  # 用于部分迹计算的累积维度
        self.state_name = state_name  # 量子态名称，用于描述该量子态

    def init_hilbert_space(self, dims: List[int], types: Optional[List[SubsystemType]] = None,
                           state_name: Optional[str] = None) -> None:
        """
        初始化希尔伯特空间的维度结构和子系统类型

        该函数设置系统的希尔伯特空间结构，记录各子系统的维度和类型（可选），并计算总维度。
        初始化后的系统尚未设置具体的量子态。

        参数：
            dims (List[int]): 子系统维度列表，例如[5, 3]表示两个子系统，维度分别为5和3
            types (List[SubsystemType], optional): 子系统类型列表，长度需与dims相同
                                                 默认为None，此时所有子系统类型为OTHER
            state_name (str, optional): 量子态名称，用于描述该量子态，默认为空字符串

        数学表示：
            假设有n个子系统，维度为 {d_1, d_2, ..., d_n}
            总维数 D_total = prod(d_j) = d_1 * d_2 * ... * d_n

        例子：
            # 初始化一个5维原子和3维腔的复合系统
            init_hilbert_space([5, 3],
                              [SubsystemType.ATOM, SubsystemType.CAVITY],
                              "原子-腔复合系统")

            # 初始化两个任意维度的子系统
            init_hilbert_space([4, 6], state_name="四能级-六能级系统")
        """
        # 验证输入维度合法性
        if not all(d > 0 and isinstance(d, int) for d in dims):
            raise ValueError("所有子系统维度必须是正整数")

        self.subsystem_dims = dims.copy()
        self.total_dim = np.prod(dims)

        # 设置子系统类型（可选）
        if types is None:
            self.subsystem_types = [SubsystemType.OTHER] * len(dims)
        else:
            if len(types) != len(dims):
                raise ValueError("子系统类型列表长度必须与维度列表长度相同")
            self.subsystem_types = types.copy()

        # 设置状态名称
        if state_name is not None:
            self.state_name = state_name
        else:
            # 自动生成状态名称
            name_parts = []
            for i, (dim, type_) in enumerate(zip(dims, self.subsystem_types)):
                if type_ == SubsystemType.ATOM:
                    name_parts.append(f"原子{i+1}({dim}态)")
                elif type_ == SubsystemType.CAVITY:
                    name_parts.append(f"腔{i+1}({dim}态)")
                else:
                    name_parts.append(f"系统{i+1}({dim}态)")
            self.state_name = "-".join(name_parts) + "复合系统"

        # 计算累积维度，用于部分迹计算
        self.subsystem_cumulative_dims = []
        cum_dim = 1
        for d in reversed(dims):
            self.subsystem_cumulative_dims.insert(0, cum_dim)
            cum_dim *= d

        # 重置量子态相关属性为None，需要显式调用set_xxx方法设置
        self.rho = None
        self.hamiltonian = None
        self.lindblad_ops = None

        # 添加基本空间信息到参数字典
        self.parameters['hilbert_dims'] = dims
        self.parameters['total_dim'] = self.total_dim
        self.parameters['subsystem_types'] = [t.value for t in self.subsystem_types]
        self.parameters['state_name'] = self.state_name

    def set_density_matrix(self, rho: Union[np.ndarray, qt.Qobj]) -> None:
        """
        设置系统的密度矩阵

        设置当前系统的密度矩阵ρ，并验证其合法性（厄米性、迹为1、半正定性）

        参数：
            rho (Union[np.ndarray, qt.Qobj]): 密度矩阵，可以是NumPy数组或QuTiP的Qobj对象

        数学表示：
            ρ_global ← ρ

        要求：
            - ρ ≥ 0 (半正定)
            - Tr(ρ) = 1 (迹为1)
            - ρ = ρ† (厄米矩阵)
            - dim(ρ) = D_total (维度匹配)

        例子：
            对于纯态 |ψ⟩，可以通过 ρ = |ψ⟩⟨ψ| 设置
            对于混合态，直接提供密度矩阵
        """
        if self.total_dim is None:
            raise ValueError("必须先调用init_hilbert_space初始化希尔伯特空间")

        # 如果是NumPy数组，转换为QuTiP的Qobj对象
        if isinstance(rho, np.ndarray):
            # 检查维度
            if rho.shape != (self.total_dim, self.total_dim):
                raise ValueError(f"密度矩阵维度应为({self.total_dim}, {self.total_dim})，"
                                f"但给定的维度为{rho.shape}")

            # 转换为QuTiP的Qobj对象
            qobj_rho = qt.Qobj(rho, dims=[self.subsystem_dims, self.subsystem_dims])
        elif isinstance(rho, qt.Qobj):
            # 检查维度
            if rho.shape != (self.total_dim, self.total_dim):
                raise ValueError(f"密度矩阵维度应为({self.total_dim}, {self.total_dim})，"
                                f"但给定的维度为{rho.shape}")

            # 确保dims参数正确
            qobj_rho = qt.Qobj(rho.full(), dims=[self.subsystem_dims, self.subsystem_dims])
        else:
            raise TypeError("密度矩阵必须是NumPy数组或QuTiP的Qobj对象")

        # 检查厄米性
        if not qobj_rho.isherm:
            raise ValueError("密度矩阵必须是厄米矩阵(ρ = ρ†)")

        # 检查迹是否等于1
        trace = qobj_rho.tr()
        if not np.isclose(trace, 1.0, atol=1e-10):
            raise ValueError(f"密度矩阵的迹必须为1，但实际值为{trace}")

        # 检查半正定性（可选，计算量较大）
        eigvals = qobj_rho.eigenenergies()
        if np.any(eigvals < -1e-10):  # 允许一些数值误差
            raise ValueError("密度矩阵必须是半正定的")

        # 存储密度矩阵
        self.rho = qobj_rho

    def get_density_matrix(self) -> qt.Qobj:
        """
        获取当前密度矩阵

        返回当前系统的密度矩阵ρ（QuTiP的Qobj对象）。

        返回：
            qt.Qobj: 当前密度矩阵（QuTiP的Qobj对象）

        数学表示：
            返回 ρ_global
        """
        if self.rho is None:
            raise ValueError(f"量子系统 '{self.state_name}' 的密度矩阵尚未初始化")
        return self.rho

    def register_parameter(self, key: str, value: Any) -> None:
        """
        注册物理参数或仿真参数

        存储与量子系统相关的物理参数或仿真参数，如波包长度τ、激发概率p(t)、衰减率γ等。

        参数：
            key (str): 参数的名称标识
            value (Any): 参数值，可以是数值、列表、字典等任何可序列化的对象

        例子：
            register_parameter('kappa', 0.5)  # 腔的衰减率κ
            register_parameter('tau', 1.0)    # 波包长度τ
            register_parameter('gamma', 0.1)  # 自发辐射率γ
        """
        self.parameters[key] = value

    def get_parameter(self, key: str) -> Any:
        """
        获取已注册的物理参数

        读取先前用register_parameter登记的参数。

        参数：
            key (str): 参数名称

        返回：
            Any: 参数值

        例子：
            gamma = get_parameter('gamma')  # 获取自发辐射率γ
        """
        if key not in self.parameters:
            raise KeyError(f"参数'{key}'未在量子系统 '{self.state_name}' 中注册")
        return self.parameters[key]

    def set_hamiltonian(self, hamiltonian: Union[np.ndarray, qt.Qobj]) -> None:
        """
        设置系统的哈密顿量

        参数：
            hamiltonian (Union[np.ndarray, qt.Qobj]): 哈密顿量矩阵，可以是NumPy数组或QuTiP的Qobj对象

        数学表示：
            H_global ← H

        要求：
            - H 是厄米矩阵 (H = H†)
            - dim(H) = D_total (维度匹配)
        """
        if self.total_dim is None:
            raise ValueError("必须先调用init_hilbert_space初始化希尔伯特空间")

        # 如果是NumPy数组，转换为QuTiP的Qobj对象
        if isinstance(hamiltonian, np.ndarray):
            # 检查维度
            if hamiltonian.shape != (self.total_dim, self.total_dim):
                raise ValueError(f"哈密顿量矩阵维度应为({self.total_dim}, {self.total_dim})，"
                                f"但给定的维度为{hamiltonian.shape}")

            # 转换为QuTiP的Qobj对象
            qobj_hamiltonian = qt.Qobj(hamiltonian, dims=[self.subsystem_dims, self.subsystem_dims])
        elif isinstance(hamiltonian, qt.Qobj):
            # 检查维度
            if hamiltonian.shape != (self.total_dim, self.total_dim):
                raise ValueError(f"哈密顿量矩阵维度应为({self.total_dim}, {self.total_dim})，"
                                f"但给定的维度为{hamiltonian.shape}")

            # 确保dims参数正确
            qobj_hamiltonian = qt.Qobj(hamiltonian.full(), dims=[self.subsystem_dims, self.subsystem_dims])
        else:
            raise TypeError("哈密顿量矩阵必须是NumPy数组或QuTiP的Qobj对象")

        # 检查厄米性
        if not qobj_hamiltonian.isherm:
            raise ValueError("哈密顿量矩阵必须是厄米矩阵(H = H†)")

        # 存储哈密顿量
        self.hamiltonian = qobj_hamiltonian

    def get_hamiltonian(self) -> qt.Qobj:
        """
        获取当前哈密顿量

        返回：
            qt.Qobj: 当前哈密顿量（QuTiP的Qobj对象）

        数学表示：
            返回 H_global
        """
        if self.hamiltonian is None:
            raise ValueError(f"量子系统 '{self.state_name}' 的哈密顿量尚未初始化")
        return self.hamiltonian

    def set_lindblad_ops(self, lindblad_ops: List[Union[np.ndarray, qt.Qobj]]) -> None:
        """
        设置系统的Lindblad算符

        参数：
            lindblad_ops (List[Union[np.ndarray, qt.Qobj]]): Lindblad算符列表，每个可以是NumPy数组或QuTiP的Qobj对象

        数学表示：
            L_global ← [L_1, L_2, ..., L_n]

        要求：
            - 每个L_i的维度为 D_total (维度匹配)
        """
        if self.total_dim is None:
            raise ValueError("必须先调用init_hilbert_space初始化希尔伯特空间")

        if not isinstance(lindblad_ops, list):
            lindblad_ops = [lindblad_ops]

        # 转换为QuTiP的Qobj对象列表
        qobj_lindblad_ops = []
        for i, op in enumerate(lindblad_ops):
            if isinstance(op, np.ndarray):
                # 检查维度
                if op.shape != (self.total_dim, self.total_dim):
                    raise ValueError(f"Lindblad算符 #{i} 维度应为({self.total_dim}, {self.total_dim})，"
                                   f"但给定的维度为{op.shape}")

                # 转换为QuTiP的Qobj对象
                qobj_op = qt.Qobj(op, dims=[self.subsystem_dims, self.subsystem_dims])
            elif isinstance(op, qt.Qobj):
                # 检查维度
                if op.shape != (self.total_dim, self.total_dim):
                    raise ValueError(f"Lindblad算符 #{i} 维度应为({self.total_dim}, {self.total_dim})，"
                                   f"但给定的维度为{op.shape}")

                # 确保dims参数正确
                qobj_op = qt.Qobj(op.full(), dims=[self.subsystem_dims, self.subsystem_dims])
            else:
                raise TypeError(f"Lindblad算符 #{i} 必须是NumPy数组或QuTiP的Qobj对象")

            qobj_lindblad_ops.append(qobj_op)

        # 存储Lindblad算符列表
        self.lindblad_ops = qobj_lindblad_ops

    def get_lindblad_ops(self) -> List[qt.Qobj]:
        """
        获取当前Lindblad算符列表

        返回：
            List[qt.Qobj]: 当前Lindblad算符列表（每个是QuTiP的Qobj对象）

        数学表示：
            返回 [L_1, L_2, ..., L_n]
        """
        if self.lindblad_ops is None:
            raise ValueError(f"量子系统 '{self.state_name}' 的Lindblad算符尚未初始化")
        return self.lindblad_ops

    def get_operator(self, op_single: qt.Qobj, target_subsystem_idx: int) -> qt.Qobj:
        """
        将单个子系统上的算符扩展到全局希尔伯特空间。

        将作用在单个子系统上的算符扩展到整个系统的张量积空间，
        在目标子系统位置放置给定算符，其他子系统位置放置单位算符。

        被调用者: NoiseModelBuilder._add_op_to_list
        调用: 无其他函数

        参数：
            op_single (qt.Qobj): 作用在单个子系统上的算符，维度必须与目标子系统维度匹配
            target_subsystem_idx (int): 目标子系统的索引

        返回：
            qt.Qobj: 扩展到全局空间的算符

        数学表示：
            对于系统 [d_1, d_2, ..., d_n]，返回
            I_1 ⊗ I_2 ⊗ ... ⊗ op_single ⊗ ... ⊗ I_n
            其中 op_single 位于索引 target_subsystem_idx 处

        注意：
            此方法功能与 core.operators.local_operator 相同，
            为了向后兼容而保留。建议新代码直接使用 local_operator 函数。
        """
        if self.subsystem_dims is None:
            raise ValueError("必须先调用init_hilbert_space初始化希尔伯特空间")

        if target_subsystem_idx < 0 or target_subsystem_idx >= len(self.subsystem_dims):
            raise ValueError(f"目标子系统索引 {target_subsystem_idx} 超出范围 [0, {len(self.subsystem_dims)-1}]")

        target_dim = self.subsystem_dims[target_subsystem_idx]
        if op_single.shape[0] != target_dim or op_single.shape[1] != target_dim:
            raise ValueError(f"算符维度 {op_single.shape} 与目标子系统维度 {target_dim} 不匹配")

        # 构建算符列表
        op_list = []
        for i, dim in enumerate(self.subsystem_dims):
            if i == target_subsystem_idx:
                op_list.append(op_single)
            else:
                op_list.append(qt.qeye(dim))

        # 返回张量积
        return qt.tensor(op_list)

    def get_state_info(self) -> Dict[str, Any]:
        """
        获取当前量子态的完整信息

        返回一个包含当前量子态所有信息的字典，便于记录系统状态帧。

        返回：
            Dict[str, Any]: 包含量子态信息的字典，包括：
                - 'name': 状态名称
                - 'dims': 子系统维度列表
                - 'types': 子系统类型列表
                - 'total_dim': 总维度
                - 'rho': 密度矩阵（NumPy数组格式，用于存储）
                - 'qutip_dims': 密度矩阵的QuTiP dims信息
                - 'hamiltonian': 哈密顿量（NumPy数组格式，用于存储）
                - 'hamiltonian_qutip_dims': 哈密顿量的QuTiP dims信息
                - 'lindblad_ops': Lindblad算符列表（NumPy数组格式，用于存储）
                - 'lindblad_qutip_dims': Lindblad算符的QuTiP dims信息
                - 'subsystem_structure': 子系统结构描述
        """
        if self.rho is None:
            raise ValueError(f"量子系统 '{self.state_name}' 的密度矩阵尚未初始化")

        # 构建子系统结构描述
        structure = []
        for i, (dim, type_) in enumerate(zip(self.subsystem_dims, self.subsystem_types)):
            if type_ == SubsystemType.ATOM:
                structure.append(f"原子{i+1}({dim}态)")
            elif type_ == SubsystemType.CAVITY:
                structure.append(f"腔{i+1}({dim}态)")
            else:
                structure.append(f"系统{i+1}({dim}态)")

        # 构建基本信息
        info = {
            'name': self.state_name,
            'dims': self.subsystem_dims.copy(),
            'types': [t.value for t in self.subsystem_types],
            'total_dim': self.total_dim,
            'rho': self.rho.full(),  # 转换为NumPy数组用于存储
            'qutip_dims': self.rho.dims,  # 保存QuTiP的dims信息
            'subsystem_structure': structure
        }

        # 添加哈密顿量信息（如果存在）
        if self.hamiltonian is not None:
            info['hamiltonian'] = self.hamiltonian.full()
            info['hamiltonian_qutip_dims'] = self.hamiltonian.dims

        # 添加Lindblad算符信息（如果存在）
        if self.lindblad_ops is not None:
            info['lindblad_ops'] = [op.full() for op in self.lindblad_ops]
            info['lindblad_qutip_dims'] = [op.dims for op in self.lindblad_ops]

        return info

# 量子态注册表，存储所有创建的量子系统
_quantum_states: Dict[str, QuantumState] = {}

def register_quantum_state(state_name: str) -> QuantumState:
    """
    注册新的量子系统

    参数：
        state_name (str): 量子系统名称

    返回：
        QuantumState: 新创建的量子系统对象
    """
    if state_name in _quantum_states:
        raise ValueError(f"量子系统名称 '{state_name}' 已存在")

    # 创建新的量子系统
    state = QuantumState(state_name)
    _quantum_states[state_name] = state

    return state

def get_quantum_state(state_name: str) -> QuantumState:
    """
    获取量子系统对象

    参数：
        state_name (str): 量子系统名称

    返回：
        QuantumState: 指定名称的量子系统对象
    """
    if state_name is None:
        raise ValueError("必须指定量子系统名称")

    # 获取指定名称的量子系统
    if state_name not in _quantum_states:
        raise ValueError(f"量子系统 '{state_name}' 不存在")

    return _quantum_states[state_name]

def list_quantum_states() -> List[str]:
    """
    列出所有已注册的量子系统名称

    返回：
        List[str]: 量子系统名称列表
    """
    return list(_quantum_states.keys())

def remove_quantum_state(state_name: str) -> None:
    """
    移除指定的量子系统

    参数：
        state_name (str): 要移除的量子系统名称
    """
    if state_name not in _quantum_states:
        raise ValueError(f"量子系统 '{state_name}' 不存在")

    # 移除量子系统
    del _quantum_states[state_name]

# 导出函数接口，使得模块调用更简洁
def init_hilbert_space(dims: List[int], types: Optional[List[SubsystemType]] = None,
                       state_name: Optional[str] = None, system_name: str = None) -> None:
    """
    初始化希尔伯特空间的维度结构和子系统类型

    包装QuantumState.init_hilbert_space方法，提供更简洁的接口。

    参数：
        dims (List[int]): 子系统维度列表
        types (List[SubsystemType], optional): 子系统类型列表（可选）
        state_name (str, optional): 量子态名称（可选）
        system_name (str): 量子系统名称

    例子：
        # 初始化一个原子-腔系统
        init_hilbert_space([5, 3],
                          [SubsystemType.ATOM, SubsystemType.CAVITY],
                          "原子-腔复合系统",
                          "system1")
    """
    quantum_state = get_quantum_state(system_name)
    quantum_state.init_hilbert_space(dims, types, state_name)

def set_density_matrix(rho: Union[np.ndarray, qt.Qobj], system_name: str = None) -> None:
    """
    设置系统的密度矩阵

    包装QuantumState.set_density_matrix方法，提供更简洁的接口。

    参数：
        rho (Union[np.ndarray, qt.Qobj]): 密度矩阵，可以是NumPy数组或QuTiP的Qobj对象
        system_name (str): 量子系统名称
    """
    quantum_state = get_quantum_state(system_name)
    quantum_state.set_density_matrix(rho)

def get_density_matrix(system_name: str = None) -> qt.Qobj:
    """
    获取当前密度矩阵

    包装QuantumState.get_density_matrix方法，提供更简洁的接口。

    参数：
        system_name (str): 量子系统名称

    返回：
        qt.Qobj: 当前密度矩阵（QuTiP的Qobj对象）
    """
    quantum_state = get_quantum_state(system_name)
    return quantum_state.get_density_matrix()

def register_parameter(key: str, value: Any, system_name: str = None) -> None:
    """
    注册物理参数或仿真参数

    包装QuantumState.register_parameter方法，提供更简洁的接口。

    参数：
        key (str): 参数的名称标识
        value (Any): 参数值
        system_name (str): 量子系统名称
    """
    quantum_state = get_quantum_state(system_name)
    quantum_state.register_parameter(key, value)

def get_parameter(key: str, system_name: str = None) -> Any:
    """
    获取已注册的物理参数

    包装QuantumState.get_parameter方法，提供更简洁的接口。

    参数：
        key (str): 参数名称
        system_name (str): 量子系统名称

    返回：
        Any: 参数值
    """
    quantum_state = get_quantum_state(system_name)
    return quantum_state.get_parameter(key)

def get_state_info(system_name: str = None) -> Dict[str, Any]:
    """
    获取当前量子态的完整信息

    包装QuantumState.get_state_info方法，提供更简洁的接口。

    参数：
        system_name (str): 量子系统名称

    返回：
        Dict[str, Any]: 包含量子态信息的字典
    """
    quantum_state = get_quantum_state(system_name)
    return quantum_state.get_state_info()

def create_frame_data(time: float, params: Dict[str, Any] = None,
                     quantum_states: List[str] = None) -> Dict[str, Any]:
    """
    创建包含多个量子系统状态信息的帧数据

    参数:
        time (float): 当前仿真时间
        params (dict, optional): 包含附加参数的字典，默认为None
        quantum_states (list, optional): 量子系统名称列表，默认为None表示包含所有系统

    返回:
        dict: 包含帧数据的字典
    """
    frame = {
        'time': time,
        'quantum_systems': []
    }

    # 确定要包含的量子系统
    if quantum_states is None:
        systems_to_include = list_quantum_states()
    else:
        systems_to_include = quantum_states

    # 添加每个量子系统的状态信息
    for system_name in systems_to_include:
        try:
            quantum_state = get_quantum_state(system_name)
            state_info = quantum_state.get_state_info()
            frame['quantum_systems'].append(state_info)
        except Exception as e:
            import warnings
            warnings.warn(f"无法获取量子系统 '{system_name}' 的状态信息: {str(e)}")

    # 添加其他附加参数
    if params:
        frame.update(params)

    return frame

# 全局函数接口
def set_hamiltonian(hamiltonian: Union[np.ndarray, qt.Qobj], system_name: str = None) -> None:
    """
    设置系统的哈密顿量

    包装QuantumState.set_hamiltonian方法，提供更简洁的接口。

    参数：
        hamiltonian (Union[np.ndarray, qt.Qobj]): 哈密顿量矩阵，可以是NumPy数组或QuTiP的Qobj对象
        system_name (str): 量子系统名称
    """
    quantum_state = get_quantum_state(system_name)
    quantum_state.set_hamiltonian(hamiltonian)

def get_hamiltonian(system_name: str = None) -> qt.Qobj:
    """
    获取当前哈密顿量

    包装QuantumState.get_hamiltonian方法，提供更简洁的接口。

    参数：
        system_name (str): 量子系统名称

    返回：
        qt.Qobj: 当前哈密顿量（QuTiP的Qobj对象）
    """
    quantum_state = get_quantum_state(system_name)
    return quantum_state.get_hamiltonian()

def set_lindblad_ops(lindblad_ops: List[Union[np.ndarray, qt.Qobj]], system_name: str = None) -> None:
    """
    设置系统的Lindblad算符

    包装QuantumState.set_lindblad_ops方法，提供更简洁的接口。

    参数：
        lindblad_ops (List[Union[np.ndarray, qt.Qobj]]): Lindblad算符列表，每个可以是NumPy数组或QuTiP的Qobj对象
        system_name (str): 量子系统名称
    """
    quantum_state = get_quantum_state(system_name)
    quantum_state.set_lindblad_ops(lindblad_ops)

def get_lindblad_ops(system_name: str = None) -> List[qt.Qobj]:
    """
    获取当前Lindblad算符列表

    包装QuantumState.get_lindblad_ops方法，提供更简洁的接口。

    参数：
        system_name (str): 量子系统名称

    返回：
        List[qt.Qobj]: 当前Lindblad算符列表（每个是QuTiP的Qobj对象）
    """
    quantum_state = get_quantum_state(system_name)
    return quantum_state.get_lindblad_ops()

if __name__ == "__main__":
    """
    自检脚本，用于测试量子系统的创建、修改和帧记录功能
    """
    import numpy as np
    import random
    import time
    try:
        from core.system_state import SystemState, record_frame, get_system_state
    except ImportError:
        from system_state import SystemState, record_frame, get_system_state

    print("===== 量子系统管理模块自检 =====")

    # 清除之前的所有量子系统和系统状态
    _quantum_states.clear()  # 清除所有已注册的量子系统
    global _system_state_instance
    _system_state_instance = None  # 重置系统状态单例

    # 创建一个新的系统状态对象，用于管理帧
    system_state = get_system_state("量子系统自检")

    # 辅助函数：打印当前存在的量子系统
    def print_current_quantum_systems():
        print(f"当前存在的量子系统: {list(_quantum_states.keys())}")

    # 步骤1：创建单原子系统（5×5密度矩阵）
    print("\n===== 1. 测试单原子系统（5×5）=====")
    single_atom = register_quantum_state("single_atom")
    init_hilbert_space([5], [SubsystemType.ATOM], "单原子系统", "single_atom")

    # 打印当前量子系统
    print_current_quantum_systems()

    # 创建初始密度矩阵
    psi_initial = qt.rand_ket(5)  # 随机纯态
    rho_initial = psi_initial * psi_initial.dag()  # 密度矩阵
    set_density_matrix(rho_initial, "single_atom")

    # 记录原始密度矩阵
    original_rho = get_density_matrix("single_atom").full()

    # 创建一个新的随机密度矩阵，而不是修改已有的
    psi_new = qt.rand_ket(5)  # 新的随机纯态
    rho_new = psi_new * psi_new.dag()  # 新的密度矩阵
    set_density_matrix(rho_new, "single_atom")
    current_rho = get_density_matrix("single_atom").full()

    # 随机选择一个元素进行比较
    i, j = random.randint(0, 4), random.randint(0, 4)

    # 验证修改
    print(f"比较位置: [{i}, {j}]")
    print(f"原始值: {original_rho[i, j]}")
    print(f"新值: {current_rho[i, j]}")
    print(f"修改是否成功: {not np.isclose(original_rho[i, j], current_rho[i, j])}")

    # 更新系统状态对象中的量子系统列表
    system_state.quantum_systems = list(_quantum_states.keys())
    print(f"系统状态中的量子系统: {system_state.quantum_systems}")

    # 记录第一帧 - 单原子系统
    print("正在记录第一帧...")
    single_atom_info = get_state_info("single_atom")  # 获取完整的状态信息
    record_frame_data = {"time": 0.1, "params": {"test": "单原子系统测试"},
                         "quantum_systems": [single_atom_info], "timestamp": datetime.datetime.now().isoformat()}
    system_state.frames.append(record_frame_data)  # 直接添加到帧列表中
    print("第一帧记录完成")

    # 步骤2：创建两个原子-腔系统（15×15密度矩阵）
    print("\n===== 2. 测试两个原子-腔系统（15×15）=====")

    # 删除先前的系统
    remove_quantum_state("single_atom")

    # 创建两个原子-腔系统
    A_cav1 = register_quantum_state("A_cav1")
    B_cav2 = register_quantum_state("B_cav2")

    # 打印当前量子系统
    print_current_quantum_systems()

    # 初始化两个原子-腔系统
    init_hilbert_space([5, 3], [SubsystemType.ATOM, SubsystemType.CAVITY], "原子A-腔1系统", "A_cav1")
    init_hilbert_space([5, 3], [SubsystemType.ATOM, SubsystemType.CAVITY], "原子B-腔2系统", "B_cav2")

    # 创建初始密度矩阵
    psi_A_initial = qt.rand_ket(15)  # 5x3=15
    rho_A_initial = psi_A_initial * psi_A_initial.dag()
    set_density_matrix(rho_A_initial, "A_cav1")

    psi_B_initial = qt.rand_ket(15)
    rho_B_initial = psi_B_initial * psi_B_initial.dag()
    set_density_matrix(rho_B_initial, "B_cav2")

    # 记录原始密度矩阵
    original_rho_A = get_density_matrix("A_cav1").full()
    original_rho_B = get_density_matrix("B_cav2").full()

    # 创建新的随机密度矩阵
    psi_A_new = qt.rand_ket(15)
    rho_A_new = psi_A_new * psi_A_new.dag()
    set_density_matrix(rho_A_new, "A_cav1")

    psi_B_new = qt.rand_ket(15)
    rho_B_new = psi_B_new * psi_B_new.dag()
    set_density_matrix(rho_B_new, "B_cav2")

    # 验证修改
    current_rho_A = get_density_matrix("A_cav1").full()
    current_rho_B = get_density_matrix("B_cav2").full()

    # 随机选择位置进行比较
    i_A, j_A = random.randint(0, 14), random.randint(0, 14)
    i_B, j_B = random.randint(0, 14), random.randint(0, 14)

    print("A_cav1 修改验证:")
    print(f"比较位置: [{i_A}, {j_A}]")
    print(f"原始值: {original_rho_A[i_A, j_A]}")
    print(f"新值: {current_rho_A[i_A, j_A]}")
    print(f"修改是否成功: {not np.isclose(original_rho_A[i_A, j_A], current_rho_A[i_A, j_A])}")

    print("\nB_cav2 修改验证:")
    print(f"比较位置: [{i_B}, {j_B}]")
    print(f"原始值: {original_rho_B[i_B, j_B]}")
    print(f"新值: {current_rho_B[i_B, j_B]}")
    print(f"修改是否成功: {not np.isclose(original_rho_B[i_B, j_B], current_rho_B[i_B, j_B])}")

    # 更新系统状态对象中的量子系统列表
    system_state.quantum_systems = list(_quantum_states.keys())
    print(f"系统状态中的量子系统: {system_state.quantum_systems}")

    # 记录第二帧 - 两个原子-腔系统
    print("正在记录第二帧...")
    A_cav1_info = get_state_info("A_cav1")  # 获取完整的状态信息
    B_cav2_info = get_state_info("B_cav2")  # 获取完整的状态信息
    record_frame_data = {"time": 0.2, "params": {"test": "两个原子-腔系统测试"},
                         "quantum_systems": [A_cav1_info, B_cav2_info], "timestamp": datetime.datetime.now().isoformat()}
    system_state.frames.append(record_frame_data)  # 直接添加到帧列表中
    print("第二帧记录完成")

    # 步骤3：创建单原子系统和双原子系统
    print("\n===== 3. 测试单原子系统（5×5）和双原子系统（25×25）=====")

    # 删除先前的系统
    remove_quantum_state("A_cav1")
    remove_quantum_state("B_cav2")

    # 创建新的系统
    single_atom = register_quantum_state("single_atom")
    double_atom = register_quantum_state("double_atom")

    # 打印当前量子系统
    print_current_quantum_systems()

    # 初始化系统
    init_hilbert_space([5], [SubsystemType.ATOM], "单原子系统", "single_atom")
    init_hilbert_space([5, 5], [SubsystemType.ATOM, SubsystemType.ATOM], "双原子系统", "double_atom")

    # 创建初始密度矩阵
    psi_single_initial = qt.rand_ket(5)
    rho_single_initial = psi_single_initial * psi_single_initial.dag()
    set_density_matrix(rho_single_initial, "single_atom")

    psi_double_initial = qt.rand_ket(25)  # 5x5=25
    rho_double_initial = psi_double_initial * psi_double_initial.dag()
    set_density_matrix(rho_double_initial, "double_atom")

    # 记录原始密度矩阵
    original_rho_single = get_density_matrix("single_atom").full()
    original_rho_double = get_density_matrix("double_atom").full()

    # 创建新的随机密度矩阵
    psi_single_new = qt.rand_ket(5)
    rho_single_new = psi_single_new * psi_single_new.dag()
    set_density_matrix(rho_single_new, "single_atom")

    psi_double_new = qt.rand_ket(25)
    rho_double_new = psi_double_new * psi_double_new.dag()
    set_density_matrix(rho_double_new, "double_atom")

    # 验证修改
    current_rho_single = get_density_matrix("single_atom").full()
    current_rho_double = get_density_matrix("double_atom").full()

    # 随机选择位置进行比较
    i_s, j_s = random.randint(0, 4), random.randint(0, 4)
    i_d, j_d = random.randint(0, 24), random.randint(0, 24)

    # 打印single_atom验证信息
    print("\nsingle_atom 修改验证:")
    print(f"比较位置: [{i_s}, {j_s}]")
    print(f"原始值: {original_rho_single[i_s, j_s]}")
    print(f"新值: {current_rho_single[i_s, j_s]}")
    print(f"修改是否成功: {not np.isclose(original_rho_single[i_s, j_s], current_rho_single[i_s, j_s])}")

    # 打印double_atom验证信息
    print("\ndouble_atom 修改验证:")
    print(f"比较位置: [{i_d}, {j_d}]")
    print(f"原始值: {original_rho_double[i_d, j_d]}")
    print(f"新值: {current_rho_double[i_d, j_d]}")
    print(f"修改是否成功: {not np.isclose(original_rho_double[i_d, j_d], current_rho_double[i_d, j_d])}")

    # 更新系统状态对象中的量子系统列表
    system_state.quantum_systems = list(_quantum_states.keys())
    print(f"系统状态中的量子系统: {system_state.quantum_systems}")

    # 记录第三帧 - 单原子系统和双原子系统
    print("正在记录第三帧...")
    single_atom_info = get_state_info("single_atom")  # 获取完整的状态信息
    double_atom_info = get_state_info("double_atom")  # 获取完整的状态信息
    record_frame_data = {"time": 0.3, "params": {"test": "单原子系统和双原子系统测试"},
                         "quantum_systems": [single_atom_info, double_atom_info], "timestamp": datetime.datetime.now().isoformat()}
    system_state.frames.append(record_frame_data)  # 直接添加到帧列表中
    print("第三帧记录完成")

    # 步骤4：创建三原子系统
    print("\n===== 4. 测试三原子系统（125×125）=====")

    # 删除先前的系统
    remove_quantum_state("single_atom")
    remove_quantum_state("double_atom")

    # 创建新的系统
    triple_atom = register_quantum_state("triple_atom")

    # 打印当前量子系统
    print_current_quantum_systems()

    # 初始化系统
    init_hilbert_space([5, 5, 5], [SubsystemType.ATOM, SubsystemType.ATOM, SubsystemType.ATOM],
                      "三原子系统", "triple_atom")

    # 创建初始随机密度矩阵 - 使用张量积来创建125x125的矩阵
    psi1_initial = qt.rand_ket(5)
    psi2_initial = qt.rand_ket(5)
    psi3_initial = qt.rand_ket(5)
    psi_triple_initial = qt.tensor(psi1_initial, psi2_initial, psi3_initial)
    rho_triple_initial = psi_triple_initial * psi_triple_initial.dag()
    set_density_matrix(rho_triple_initial, "triple_atom")

    # 记录原始密度矩阵的特定元素
    original_rho_triple = get_density_matrix("triple_atom")

    # 使用特定位置的元素测试，而不是生成完整的密度矩阵数组
    test_i, test_j = random.randint(0, 124), random.randint(0, 124)
    original_value = original_rho_triple[test_i, test_j]

    # 创建一个新的随机态，而不是直接修改矩阵元素
    psi1_new = qt.rand_ket(5)
    psi2_new = qt.rand_ket(5)
    psi3_new = qt.rand_ket(5)
    psi_triple_new = qt.tensor(psi1_new, psi2_new, psi3_new)
    rho_triple_new = psi_triple_new * psi_triple_new.dag()
    set_density_matrix(rho_triple_new, "triple_atom")

    # 验证修改
    current_rho_triple = get_density_matrix("triple_atom")
    current_value = current_rho_triple[test_i, test_j]

    print(f"测试位置: [{test_i}, {test_j}]")
    print(f"原始值: {original_value}")
    print(f"修改后: {current_value}")
    print(f"修改是否成功: {not np.isclose(original_value, current_value)}")

    # 更新系统状态对象中的量子系统列表
    system_state.quantum_systems = list(_quantum_states.keys())
    print(f"系统状态中的量子系统: {system_state.quantum_systems}")

    # 记录第四帧 - 三原子系统
    print("正在记录第四帧...")
    triple_atom_info = get_state_info("triple_atom")  # 获取完整的状态信息
    record_frame_data = {"time": 0.4, "params": {"test": "三原子系统测试"},
                         "quantum_systems": [triple_atom_info], "timestamp": datetime.datetime.now().isoformat()}
    system_state.frames.append(record_frame_data)  # 直接添加到帧列表中
    print("第四帧记录完成")

    # 在帧信息摘要部分之前睡眠一小段时间，确保所有帧都已记录
    time.sleep(0.1)

    # 打印帧信息摘要之前先进行一些分隔
    print("\n" + "="*50)

    # 帧信息摘要
    print("\n===== 帧信息摘要 =====")
    frames = system_state.frames
    print(f"共有 {len(frames)} 个帧")

    # 设置NumPy打印选项，控制显示精度为小数点后两位
    np.set_printoptions(precision=2, suppress=True)

    for i, frame in enumerate(frames, 1):
        print(f"\n帧 {i}:")
        print(f"帧键: {list(frame.keys())}")
        if 'time' in frame:
            print(f"仿真时间: {frame['time']}")

        # 参数信息显示
        if 'params' in frame and isinstance(frame['params'], dict):
            params = frame['params']
            print(f"参数键: {list(params.keys())}")
            for key, value in params.items():
                print(f"{key}: {value}")

        # 量子系统信息显示
        if 'quantum_systems' in frame and isinstance(frame['quantum_systems'], list):
            quantum_systems = frame['quantum_systems']
            print(f"  共有 {len(quantum_systems)} 个量子系统")

            # 遍历每个量子系统并显示其信息
            for idx, sys_info in enumerate(quantum_systems):
                if not isinstance(sys_info, dict):
                    continue

                print(f"  - 量子系统 {idx+1}:")
                if 'name' in sys_info:
                    print(f"    名称: {sys_info['name']}")
                if 'total_dim' in sys_info:
                    print(f"    总维度: {sys_info['total_dim']}")
                if 'dims' in sys_info:
                    print(f"    子系统维度: {sys_info['dims']}")
                if 'types' in sys_info:
                    print(f"    子系统类型: {sys_info['types']}")
                if 'subsystem_structure' in sys_info:
                    print(f"    子系统结构: {sys_info['subsystem_structure']}")

                # 打印密度矩阵（仅显示到小数点后两位）
                if 'rho' in sys_info and isinstance(sys_info['rho'], np.ndarray):
                    rho = sys_info['rho']
                    print(f"    密度矩阵 (维度 {rho.shape}):")
                    # 对于大型矩阵，只打印左上角10x10的子矩阵
                    if rho.shape[0] > 10:
                        print("    (仅显示左上角10x10子矩阵)")
                        print(rho[:10, :10])
                    else:
                        print(rho)

        print("-" * 40)

    # 恢复NumPy默认打印选项
    np.set_printoptions(precision=8, suppress=False)

    print("\n===== 自检完成 =====")
