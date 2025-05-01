# coding=utf-8
# operators.py

"""
量子算符构建模块 (基于 QuTiP)

本模块提供灵活的函数来构建和操作量子算符，特别适用于多能级原子和腔系统。
通过这些函数，可以方便地构建复杂的哈密顿量和Lindblad算符。

主要功能：
1. 将小维度算符嵌入到更大的希尔伯特空间中
2. 将单体算符扩展到多体系统的张量积空间
3. 构造多体相互作用项
4. 创建原子系统的常用算符、哈密顿量和Lindblad项
5. 创建腔场系统的常用算符、哈密顿量和Lindblad项
6. 创建原子-腔相互作用的算符
7. 辅助函数：构建总哈密顿量、时间依赖哈密顿量、旋转坐标系变换等

函数调用关系图：
```mermaid
graph TD
    %% 基础量子算符构建函数
    embed_op[embed_op_in_subspace] --> None
    local_op[local_operator] --> None
    interact[interaction_term] --> local_op

    %% 集成量子算符生成函数
    create_atom[create_atomic_operators] --> embed_op
    create_atom --> interact
    create_cavity[create_cavity_operators] --> None
    create_interact[create_interaction_operators] --> None

    %% 辅助函数
    build_total[build_total_hamiltonian] --> None
    build_time[build_time_dependent_hamiltonian] --> None
    apply_rotate[apply_rotating_frame] --> None
    check_herm[check_hermitian] --> None
```
"""

import numpy as np
import qutip as qt
from typing import List, Dict, Optional, Union, Callable, Any, Tuple


# ================ 量子算符构建函数 ================

def embed_op_in_subspace(op_small: qt.Qobj,
                         index_map: List[int],
                         big_dim: int) -> qt.Qobj:
    """
    将仅作用于少数能级的算符 op_small（形状 m×m）
    嵌入到一个 big_dim×big_dim 的更大希尔伯特空间中。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    op_small : qutip.Qobj
        小算符 (m×m)，如针对两能级的 [[0,Ω/2],[Ω/2, -Δ]]。
    index_map : list of int
        指明旧算符 op_small 的第 row/col 索引 (0..m-1)
        应映射到新空间的哪几个能级 (0..big_dim-1)。
        例如如果 op_small 是 2×2，index_map = [0,2]
        就表示旧 |0>,|1> 分别对应新空间的 |0>,|2>。
    big_dim : int
        新空间的维度（如 5 表示五能级原子）。

    返回:
    -------
    Qobj, 形状 (big_dim×big_dim),
    其在 index_map 对应的子空间里与 op_small 完全相同，
    其他行列元素全为 0。

    示例:
    -------
    >>> Delta = 1.0
    >>> Omega = 0.2
    >>> H_2lvl = qt.Qobj([[-Delta, Omega/2], [Omega/2, 0]])
    >>> embed_idx_map = [0, 2]  # 旧 |g>=0 -> 新系统的 0级; 旧 |e>=1 -> 新系统的 2级
    >>> H_5lvl_drive = embed_op_in_subspace(H_2lvl, embed_idx_map, 5)
    """
    # 先用 NumPy 矩阵来装（也可使用稀疏格式）
    data = np.zeros((big_dim, big_dim), dtype=complex)

    # 嵌入对应的子块
    for r_old, r_new in enumerate(index_map):
        for c_old, c_new in enumerate(index_map):
            data[r_new, c_new] = op_small[r_old, c_old]

    # 转成 Qobj，维度标签为 [[big_dim], [big_dim]]
    return qt.Qobj(data, dims=[[big_dim], [big_dim]])


def local_operator(op_single: qt.Qobj,
                   who: int,
                   dims: List[int]) -> qt.Qobj:
    """
    将单个子系统上的算符 op_single 嵌入到多子系统直积空间。

    被调用者: 外部调用, interaction_term
    调用: 无其他函数

    参数:
    -------
    op_single : qutip.Qobj
        作用在单个子系统上的算符，维度必须与 dims[who] 相匹配。
    who : int
        指定这是作用在第几个子系统 (0-based) 上的算符。
    dims : list of int
        整个系统的各子空间维度, 如 [5, 3, 5] 表示三部分:
         - 第0个子系统5维(五能级原子)
         - 第1个子系统3维(单模腔)
         - 第2个子系统5维(另一个五能级原子)

    返回:
    -------
    Qobj, 其尺寸是 (Πdims × Πdims)，
    在第 who 个子空间放置 op_single，其他子空间放置单位阵。

    示例:
    -------
    >>> H_5lvl_drive_single = embed_op_in_subspace(H_2lvl, [0, 2], 5)
    >>> H_drive_atom1 = local_operator(H_5lvl_drive_single, who=1, dims=[5,5,5])
    """
    op_list = []
    for i, d in enumerate(dims):
        if i == who:
            op_list.append(op_single)
        else:
            op_list.append(qt.qeye(d))
    return qt.tensor(op_list)


def interaction_term(op_i: qt.Qobj,
                     i: int,
                     op_j: qt.Qobj,
                     j: int,
                     factor: float,
                     dims: List[int]) -> qt.Qobj:
    """
    构造两体相互作用项 (如里德堡阻塞、偶极-偶极等)，
    只在第 i 和第 j 个子系统上非平凡，其他子系统为单位算符。

    被调用者: 外部调用
    调用: local_operator

    参数:
    -------
    op_i : qutip.Qobj
        作用在第 i 个子系统上的算符(维度 dims[i]×dims[i])。
    i : int
        子系统 i 的索引。
    op_j : qutip.Qobj
        作用在第 j 个子系统上的算符(维度 dims[j]×dims[j])。
    j : int
        子系统 j 的索引。
    factor : float
        该相互作用的整体系数 (coupling strength, V_rr 等)。
    dims : list of int
        整个系统的各子空间维度。

    返回:
    -------
    Qobj, 其尺寸是 (Πdims × Πdims)。
    表达式形如 factor * (I⊗...⊗op_i⊗...⊗op_j⊗...)。

    示例:
    -------
    >>> P_r = qt.basis(5,4)*qt.basis(5,4).dag()  # 5x5, 在 |4> 上投影
    >>> V_rr = 2*np.pi*10e-3  # 10 MHz
    >>> H_block_01 = interaction_term(P_r, 0, P_r, 1, V_rr, dims=[5,5,5])
    """
    assert i != j, "interaction_term: i 与 j 必须是不同子系统！"
    op_list = []
    for idx, d in enumerate(dims):
        if   idx == i: op_list.append(op_i)
        elif idx == j: op_list.append(op_j)
        else:          op_list.append(qt.qeye(d))
    return factor * qt.tensor(op_list)


# ================ 集成量子算符生成函数 ================

def create_atomic_operators(dim: int, level_labels: Optional[List[str]] = None,
                           energies: Optional[List[float]] = None) -> Dict[str, Any]:
    """
    创建原子系统的常用算符，包括投影算符、跃迁算符、哈密顿量和Lindblad项。

    被调用者: 外部调用
    调用: embed_op_in_subspace, interaction_term

    参数:
    -------
    dim : int
        原子系统的希尔伯特空间维度
    level_labels : Optional[List[str]]
        能级标签列表，如 ['g', 's', 'e0', 'e1', 'r']
    energies : Optional[List[float]]
        各能级的能量，长度应等于dim。若为None，则所有能级能量为0

    返回:
    -------
    Dict[str, Any]: 常用算符的字典，包括：
        - 基本算符:
          - 'P_i': 第i个能级的投影算符 |i⟩⟨i|
          - 'sigma_ij': 从j到i的跃迁算符 |i⟩⟨j|
          - 对于二能级系统，还包括泡利算符 'sigma_x', 'sigma_y', 'sigma_z'
        - 哈密顿量:
          - 'H_0': 本征能量哈密顿量
        - 构建函数:
          - 'drive': 函数，用于构建驱动哈密顿量
          - 'spontaneous_emission': 函数，用于构建自发辐射Lindblad算符
          - 'dephasing': 函数，用于构建退相位Lindblad算符
          - 'rydberg_interaction': 函数，用于构建Rydberg相互作用哈密顿量
    """
    ops = {}

    # 使用默认标签或提供的标签
    if level_labels is None:
        if dim == 2:
            level_labels = ['g', 'e']
        elif dim == 3:
            level_labels = ['g', 's', 'e']
        elif dim == 4:
            level_labels = ['g', 's', 'e0', 'e1']
        elif dim == 5:
            level_labels = ['g', 's', 'e0', 'e1', 'r']
        else:
            level_labels = [str(i) for i in range(dim)]

    # 投影算符
    for i in range(dim):
        # 数字索引的投影算符
        ops[f'P_{i}'] = qt.basis(dim, i) * qt.basis(dim, i).dag()
        # 使用能级标签的投影算符
        if i < len(level_labels):
            ops[f'P_{level_labels[i]}'] = ops[f'P_{i}']

    # 跃迁算符
    for i in range(dim):
        for j in range(dim):
            if i != j:
                # 数字索引的跃迁算符
                ops[f'sigma_{i}{j}'] = qt.basis(dim, i) * qt.basis(dim, j).dag()
                # 使用能级标签的跃迁算符
                if i < len(level_labels) and j < len(level_labels):
                    ops[f'sigma_{level_labels[i]}{level_labels[j]}'] = ops[f'sigma_{i}{j}']

    # 对于二能级系统，添加泡利算符
    if dim == 2:
        ops['sigma_x'] = qt.sigmax()
        ops['sigma_y'] = qt.sigmay()
        ops['sigma_z'] = qt.sigmaz()
        ops['sigma_plus'] = qt.sigmap()
        ops['sigma_minus'] = qt.sigmam()

    # 本征能量哈密顿量
    if energies is None:
        energies = np.zeros(dim)
    elif len(energies) != dim:
        raise ValueError(f"能量列表长度 {len(energies)} 与系统维度 {dim} 不匹配")

    ops['H_0'] = qt.Qobj(np.diag(energies), dims=[[dim], [dim]])

    # 驱动哈密顿量构建函数
    def drive(rabi_freq: float, detuning: float, level_i: int, level_j: int) -> qt.Qobj:
        """
        构建描述经典驱动场（激光/微波）的哈密顿量。

        参数:
        -------
        rabi_freq : float
            Rabi频率（单位：rad/s）
        detuning : float
            失谐（单位：rad/s）
        level_i : int
            第一个能级索引
        level_j : int
            第二个能级索引（通常是更高能级）

        返回:
        -------
        qt.Qobj: 驱动哈密顿量，形式为 (Ω/2)(|i⟩⟨j|+|j⟩⟨i|) - Δ|j⟩⟨j|
        """
        # 添加耦合项
        H_coupling = (rabi_freq/2) * (ops[f'sigma_{level_i}{level_j}'] + ops[f'sigma_{level_j}{level_i}'])

        # 添加失谐项
        H_detuning = -detuning * ops[f'P_{level_j}']

        return H_coupling + H_detuning

    ops['drive'] = drive

    # 自发辐射Lindblad算符构建函数
    def spontaneous_emission(gamma: float, level_g: int, level_e: int) -> qt.Qobj:
        """
        构建描述原子自发辐射的Lindblad算符。

        参数:
        -------
        gamma : float
            自发辐射速率（单位：rad/s）
        level_g : int
            基态索引
        level_e : int
            激发态索引

        返回:
        -------
        qt.Qobj: 自发辐射算符，形式为 √γ|g⟩⟨e|
        """
        return np.sqrt(gamma) * ops[f'sigma_{level_g}{level_e}']

    ops['spontaneous_emission'] = spontaneous_emission

    # 退相位Lindblad算符构建函数
    def dephasing(gamma_phi: float, level: int) -> qt.Qobj:
        """
        构建描述纯退相位的Lindblad算符。

        参数:
        -------
        gamma_phi : float
            退相位速率（单位：rad/s）
        level : int
            受影响能级的索引

        返回:
        -------
        qt.Qobj: 退相位算符，形式为 √γ_φ|level⟩⟨level|
        """
        return np.sqrt(gamma_phi) * ops[f'P_{level}']

    ops['dephasing'] = dephasing

    # Rydberg相互作用哈密顿量构建函数
    def rydberg_interaction(V_rr: float, r_level: int, dims: List[int]) -> qt.Qobj:
        """
        构建描述两个原子之间Rydberg相互作用的哈密顿量。

        参数:
        -------
        V_rr : float
            Rydberg相互作用强度（单位：rad/s）
        r_level : int
            Rydberg态的索引
        dims : List[int]
            系统的维度列表，如 [dim, dim] 表示两个相同的原子

        返回:
        -------
        qt.Qobj: Rydberg相互作用哈密顿量，形式为 V_rr|r⟩₁⟨r|⊗|r⟩₂⟨r|
        """
        # 使用interaction_term函数构建相互作用项
        return interaction_term(ops[f'P_{r_level}'], 0, ops[f'P_{r_level}'], 1, V_rr, dims)

    ops['rydberg_interaction'] = rydberg_interaction

    return ops


def create_cavity_operators(dim: int, cavity_freq: Optional[float] = None) -> Dict[str, Any]:
    """
    创建腔场系统的常用算符，包括基本算符、哈密顿量和Lindblad项。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    dim : int
        腔场的截断维度（最大光子数+1）
    cavity_freq : Optional[float]
        腔场频率（单位：rad/s），若为None则设为0

    返回:
    -------
    Dict[str, Any]: 常用算符的字典，包括：
        - 基本算符:
          - 'a': 湮灭算符
          - 'a_dag': 创生算符
          - 'n': 数算符
          - 'P_i': 光子数态投影算符 |i⟩⟨i|
        - 哈密顿量:
          - 'H_0': 腔场能量哈密顿量 ħω_c a†a
        - 构建函数:
          - 'decay': 函数，用于构建腔场衰减Lindblad算符
          - 'displacement': 函数，用于构建相干驱动哈密顿量
    """
    ops = {}

    # 基本算符
    ops['a'] = qt.destroy(dim)
    ops['a_dag'] = qt.create(dim)
    ops['n'] = ops['a_dag'] * ops['a']

    # 光子数态投影
    for i in range(dim):
        ops[f'P_{i}'] = qt.basis(dim, i) * qt.basis(dim, i).dag()

    # 腔场能量哈密顿量
    if cavity_freq is None:
        cavity_freq = 0.0

    ops['H_0'] = cavity_freq * ops['n']

    # 腔场衰减Lindblad算符构建函数
    def decay(kappa: float) -> qt.Qobj:
        """
        构建描述腔场衰减的Lindblad算符。

        参数:
        -------
        kappa : float
            腔场衰减速率（单位：rad/s）

        返回:
        -------
        qt.Qobj: 腔场衰减算符，形式为 √κ·a
        """
        return np.sqrt(kappa) * ops['a']

    ops['decay'] = decay

    # 腔场相干驱动哈密顿量构建函数
    def displacement(epsilon: float, phase: float = 0.0) -> qt.Qobj:
        """
        构建描述腔场相干驱动的哈密顿量。

        参数:
        -------
        epsilon : float
            驱动强度（单位：rad/s）
        phase : float
            驱动相位（单位：rad）

        返回:
        -------
        qt.Qobj: 驱动哈密顿量，形式为 ε(e^{iφ}a† + e^{-iφ}a)
        """
        return epsilon * (np.exp(1j*phase) * ops['a_dag'] + np.exp(-1j*phase) * ops['a'])

    ops['displacement'] = displacement

    return ops


def create_interaction_operators(atom_ops: Dict[str, Any], cavity_ops: Dict[str, Any],
                               cavity_dim: int) -> Dict[str, Any]:
    """
    创建描述原子-腔相互作用的算符。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    atom_ops : Dict[str, Any]
        原子系统算符字典，由create_atomic_operators生成
    cavity_ops : Dict[str, Any]
        腔场系统算符字典，由create_cavity_operators生成
    cavity_dim : int
        腔场的截断维度

    返回:
    -------
    Dict[str, Any]: 相互作用算符的字典，包括：
        - 构建函数:
          - 'jaynes_cummings': 函数，用于构建Jaynes-Cummings哈密顿量
          - 'beam_splitter': 函数，用于构建分束器哈密顿量
    """
    ops = {}

    # Jaynes-Cummings哈密顿量构建函数
    def jaynes_cummings(g_coupling: float, level_g: int, level_e: int) -> qt.Qobj:
        """
        构建Jaynes-Cummings哈密顿量，描述原子-腔耦合。

        参数:
        -------
        g_coupling : float
            原子-腔耦合强度（单位：rad/s）
        level_g : int
            基态索引
        level_e : int
            激发态索引

        返回:
        -------
        qt.Qobj: Jaynes-Cummings哈密顿量，形式为 g(σ⁺a + σ⁻a†)
        """
        # 获取原子跃迁算符
        sigma_plus = atom_ops[f'sigma_{level_e}{level_g}']
        sigma_minus = atom_ops[f'sigma_{level_g}{level_e}']

        # 构建相互作用项
        H_int = g_coupling * (qt.tensor(sigma_plus, cavity_ops['a']) +
                             qt.tensor(sigma_minus, cavity_ops['a_dag']))

        return H_int

    ops['jaynes_cummings'] = jaynes_cummings

    # 分束器哈密顿量构建函数
    def beam_splitter(kappa: float) -> qt.Qobj:
        """
        构建描述两个腔模式之间分束器耦合的哈密顿量。

        参数:
        -------
        kappa : float
            分束器耦合强度（单位：rad/s）

        返回:
        -------
        qt.Qobj: 分束器哈密顿量，形式为 κ(a₁†a₂ + a₂†a₁)
        """
        # 创建两个腔场的算符
        a1 = qt.tensor(cavity_ops['a'], qt.qeye(cavity_dim))
        a2 = qt.tensor(qt.qeye(cavity_dim), cavity_ops['a'])

        # 构建分束器哈密顿量
        H_BS = kappa * (a1.dag() * a2 + a2.dag() * a1)

        return H_BS

    ops['beam_splitter'] = beam_splitter

    return ops


# ================ 辅助函数 ================

def build_total_hamiltonian(h_terms: List[qt.Qobj]) -> qt.Qobj:
    """
    将多个哈密顿量项组合成总哈密顿量。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    h_terms : List[qt.Qobj]
        哈密顿量项列表

    返回:
    -------
    qt.Qobj: 总哈密顿量
    """
    if not h_terms:
        raise ValueError("哈密顿量项列表不能为空")

    # 检查所有项的维度是否一致
    dims = h_terms[0].dims
    for i, h in enumerate(h_terms):
        if h.dims != dims:
            raise ValueError(f"哈密顿量项 {i} 的维度 {h.dims} 与第一项的维度 {dims} 不匹配")

    # 相加所有项
    h_total = sum(h_terms)
    return h_total


def build_time_dependent_hamiltonian(h_list: List[qt.Qobj],
                                    coeffs: List[Callable[[float], float]]) -> List:
    """
    构建时间依赖哈密顿量，适用于QuTiP的求解器。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    h_list : List[qt.Qobj]
        哈密顿量项列表
    coeffs : List[Callable[[float], float]]
        时间依赖系数函数列表，与h_list一一对应

    返回:
    -------
    List: QuTiP格式的时间依赖哈密顿量 [h_list, coeffs]
    """
    if len(h_list) != len(coeffs):
        raise ValueError(f"哈密顿量项数量 {len(h_list)} 与系数函数数量 {len(coeffs)} 不匹配")

    return [h_list, coeffs]


def apply_rotating_frame(hamiltonian: qt.Qobj,
                         frame_operator: qt.Qobj,
                         omega: float) -> qt.Qobj:
    """
    将哈密顿量变换到旋转坐标系。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    hamiltonian : qt.Qobj
        原始哈密顿量
    frame_operator : qt.Qobj
        定义旋转坐标系的算符，如 σz/2 或 a†a
    omega : float
        旋转频率（单位：rad/s）

    返回:
    -------
    qt.Qobj: 旋转坐标系下的哈密顿量 H_rot = H - ω·A
    """
    return hamiltonian - omega * frame_operator


def check_hermitian(operator: qt.Qobj, tol: float = 1e-10) -> bool:
    """
    检查算符是否是厄米的（自伴的）。

    被调用者: 外部调用
    调用: 无其他函数

    参数:
    -------
    operator : qt.Qobj
        要检查的算符
    tol : float
        容差

    返回:
    -------
    bool: 如果算符是厄米的，则为True
    """
    diff = operator - operator.dag()
    return diff.norm() < tol


# ================ 测试代码 ================

def main():
    """
    测试量子算符构建函数的主函数。

    分三步测试：
    1. 构建5能级原子的哈密顿量和Lindblad项
    2. 构建3能级光学腔的哈密顿量
    3. 将原子和腔系统组合成15×15的哈密顿量并添加耦合项
    """
    print("=" * 80)
    print("测试量子算符构建函数")
    print("=" * 80)

    # 步骤1：构建5能级原子的哈密顿量和Lindblad项
    print("\n步骤1：构建5能级原子的哈密顿量和Lindblad项")
    print("-" * 60)

    # 创建5能级原子系统
    atom_dim = 5
    level_labels = ['g', 's', 'e0', 'e1', 'r']
    # 能级能量（单位：GHz，相对于基态）
    energies = [0.0, 1.0, 2.0, 2.1, 10.0]

    print(f"创建{atom_dim}能级原子系统，能级标签：{level_labels}")
    print(f"能级能量（GHz）：{energies}")

    # 创建原子算符
    atom_ops = create_atomic_operators(atom_dim, level_labels, energies)

    # 构建原子哈密顿量
    H_atom = atom_ops['H_0']
    print("\n原子本征能量哈密顿量：")
    print(H_atom)

    # 添加驱动项（|g⟩ ↔ |e0⟩）
    rabi_freq = 0.1  # Rabi频率（GHz）
    detuning = 0.05  # 失谐（GHz）
    H_drive = atom_ops['drive'](rabi_freq, detuning, 0, 2)  # |g⟩(0) ↔ |e0⟩(2)
    print(f"\n添加驱动项 |g⟩ ↔ |e0⟩，Rabi频率：{rabi_freq} GHz，失谐：{detuning} GHz")
    print(H_drive)

    # 构建总原子哈密顿量
    H_atom_total = H_atom + H_drive
    print("\n总原子哈密顿量：")
    print(H_atom_total)

    # 构建Lindblad项
    gamma = 0.01  # 自发辐射速率（GHz）
    L_spontaneous = atom_ops['spontaneous_emission'](gamma, 0, 2)  # |e0⟩(2) → |g⟩(0)
    print(f"\n自发辐射Lindblad项 |e0⟩ → |g⟩，速率：{gamma} GHz")
    print(L_spontaneous)

    # 添加退相位
    gamma_phi = 0.005  # 退相位速率（GHz）
    L_dephasing = atom_ops['dephasing'](gamma_phi, 4)  # Rydberg态退相位
    print(f"\nRydberg态退相位Lindblad项，速率：{gamma_phi} GHz")
    print(L_dephasing)

    # 验证哈密顿量是否是厄米的
    is_hermitian = check_hermitian(H_atom_total)
    print(f"\n总原子哈密顿量是否厄米：{is_hermitian}")

    # 步骤2：构建3能级光学腔的哈密顿量
    print("\n\n步骤2：构建3能级光学腔的哈密顿量")
    print("-" * 60)

    # 创建3能级光学腔系统（0, 1, 2光子态）
    cavity_dim = 3
    cavity_freq = 2.0  # 腔频率（GHz）

    print(f"创建{cavity_dim}能级光学腔系统（0, 1, 2光子态）")
    print(f"腔频率：{cavity_freq} GHz")

    # 创建腔场算符
    cavity_ops = create_cavity_operators(cavity_dim, cavity_freq)

    # 构建腔场哈密顿量
    H_cavity = cavity_ops['H_0']
    print("\n腔场能量哈密顿量：")
    print(H_cavity)

    # 添加相干驱动
    epsilon = 0.05  # 驱动强度（GHz）
    H_displacement = cavity_ops['displacement'](epsilon)
    print(f"\n添加相干驱动，强度：{epsilon} GHz")
    print(H_displacement)

    # 构建总腔场哈密顿量
    H_cavity_total = H_cavity + H_displacement
    print("\n总腔场哈密顿量：")
    print(H_cavity_total)

    # 构建腔场衰减Lindblad项
    kappa = 0.02  # 腔衰减速率（GHz）
    L_cavity = cavity_ops['decay'](kappa)
    print(f"\n腔场衰减Lindblad项，速率：{kappa} GHz")
    print(L_cavity)

    # 验证哈密顿量是否是厄米的
    is_hermitian = check_hermitian(H_cavity_total)
    print(f"\n总腔场哈密顿量是否厄米：{is_hermitian}")

    # 步骤3：将原子和腔系统组合成15×15的哈密顿量并添加耦合项
    print("\n\n步骤3：将原子和腔系统组合成15×15的哈密顿量并添加耦合项")
    print("-" * 60)

    print(f"组合{atom_dim}能级原子和{cavity_dim}能级腔场系统")
    print(f"总维度：{atom_dim} × {cavity_dim} = {atom_dim * cavity_dim}")

    # 创建相互作用算符
    interaction_ops = create_interaction_operators(atom_ops, cavity_ops, cavity_dim)

    # 构建复合系统哈密顿量
    # 1. 原子部分
    H_atom_in_composite = qt.tensor(H_atom_total, qt.qeye(cavity_dim))
    print("\n复合系统中的原子哈密顿量：")
    print(f"维度：{H_atom_in_composite.shape}")

    # 2. 腔场部分
    H_cavity_in_composite = qt.tensor(qt.qeye(atom_dim), H_cavity_total)
    print("\n复合系统中的腔场哈密顿量：")
    print(f"维度：{H_cavity_in_composite.shape}")

    # 3. Jaynes-Cummings相互作用
    g_coupling = 0.1  # 原子-腔耦合强度（GHz）
    H_JC = interaction_ops['jaynes_cummings'](g_coupling, 0, 2)  # |g⟩(0) ↔ |e0⟩(2)
    print(f"\nJaynes-Cummings相互作用，耦合强度：{g_coupling} GHz")
    print(f"维度：{H_JC.shape}")

    # 构建总哈密顿量
    H_total = H_atom_in_composite + H_cavity_in_composite + H_JC
    print("\n总系统哈密顿量：")
    print(f"维度：{H_total.shape}")
    print("哈密顿量内容：")
    print(H_total)

    # 构建复合系统的Lindblad项
    # 1. 原子自发辐射
    L_spontaneous_in_composite = qt.tensor(L_spontaneous, qt.qeye(cavity_dim))
    print("\n复合系统中的原子自发辐射Lindblad项：")
    print(f"维度：{L_spontaneous_in_composite.shape}")

    # 2. 腔场衰减
    L_cavity_in_composite = qt.tensor(qt.qeye(atom_dim), L_cavity)
    print("\n复合系统中的腔场衰减Lindblad项：")
    print(f"维度：{L_cavity_in_composite.shape}")

    # 验证总哈密顿量是否是厄米的
    is_hermitian = check_hermitian(H_total)
    print(f"\n总系统哈密顿量是否厄米：{is_hermitian}")

    print("\n测试完成！")


if __name__ == "__main__":
    main()