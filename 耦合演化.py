import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from typing import List, Dict
import scipy.sparse as sp

# === 第一部分：TwoCavityThreeLevelAtom ===
class TwoCavityThreeLevelAtom:
    """
    两腔耦合的三能级原子-QED模型：
      - 每个光学腔光子数0/1（fock_dim=2）
      - 原子3能级
      - 子系统Hilbert空间=2×3=6
      - Hamiltonian 显式使用 λ, L, w0, F, FSR, Δν, V, g_max, κ 等物理量
      - 构造有效非厄米 Hamiltonian:
        H_eff = ħ [ω_c a†a ⊗ I + I ⊗ H_atom + g_max (a†⊗σ10 + a⊗σ10†)]
              − i ħ (κ/2) a†a ⊗ I
      - 考虑腔-腔耦合，通过 Schrieffer–Wolff 变换得到 5×5 原子有效哈密顿量
    """
    c = 2.99792458e8      
    ħ = 1.0545718e-34     

    def __init__(self,
                 cavity_configs: List[Dict[str, float]],
                 atom_levels: List[float],
                 cavity_coupling: float = 0.0):
        assert len(cavity_configs) == 2, "需要两个腔的配置"
        assert len(atom_levels) == 3,   "需要三能级原子"
        self.atom_levels = atom_levels  
        self.g_cav = cavity_coupling     

        self.d_cav  = 2
        self.d_atom = 3
        self.a       = qt.destroy(self.d_cav)
        self.I_cav   = qt.qeye(self.d_cav)
        self.I_atom  = qt.qeye(self.d_atom)

        self.proj = [qt.basis(self.d_atom, i) * qt.basis(self.d_atom, i).dag()
                     for i in range(self.d_atom)]
        self.sigma10 = qt.basis(self.d_atom, 0) * qt.basis(self.d_atom, 1).dag()

        self.params = []
        for cfg in cavity_configs:
            lam   = cfg["wavelength"]       
            L     = cfg["cavity_length"]    
            w0    = cfg["cavity_waist"]     
            F     = cfg["finesse"]          
            γ_sp  = cfg["gamma_sp"]         

            fsr   = self.c / (2 * L)          
            Δν    = fsr / F                   
            Vmod  = np.pi * w0**2 * L / 4     
            g_max = np.sqrt(3 * lam**2 * γ_sp / (8 * np.pi * Vmod))
            ω_c   = 2 * np.pi * (self.c / lam) 
            κ     = np.pi * Δν                

            self.params.append({
                "wavelength": lam,
                "cavity_length": L,
                "cavity_waist": w0,
                "finesse": F,
                "FSR": fsr,
                "linewidth": Δν,
                "mode_volume": Vmod,
                "g_max": g_max,
                "ω_c": ω_c,
                "κ": κ
            })

    def generate_effective_hamiltonians(self) -> List[qt.Qobj]:
        """
        返回两个 6×6 非厄米 H_eff，用于单腔独立演化
        """
        E0, E1, E2 = self.atom_levels
        H_atom = (E0 * self.proj[0] + E1 * self.proj[1] + E2 * self.proj[2])

        H_list = []
        for p in self.params:
            H_cav   = p["ω_c"] * qt.tensor(self.a.dag() * self.a, self.I_atom)
            H_at    = qt.tensor(self.I_cav, H_atom)
            H_int   = p["g_max"] * (
                          qt.tensor(self.a.dag(), self.sigma10)
                        + qt.tensor(self.a, self.sigma10.dag())
                      )
            H_decay = -1j * (p["κ"] / 2) * qt.tensor(self.a.dag() * self.a,
                                                     self.I_atom)
            H_eff = self.ħ * (H_cav + H_at + H_int) + H_decay
            H_list.append(H_eff)
        return H_list

    def generate_sw_effective_hamiltonian(self) -> qt.Qobj:
        """
        Schrieffer–Wolff 变换：在零光子子空间构造 5×5 原子有效哈密顿量
          - 基底顺序：|0,0>, |1,0>, |0,1>, |2,0>, |0,2>
          - 包含单体拉姆位移 δ_n、原子间交换耦合 j 及腔-腔耦合 g_cav
        """
        E0, E1, E2 = self.atom_levels
        g_vals  = [p['g_max'] for p in self.params]
        wc_vals = [p['ω_c']   for p in self.params]
        omega   = [(E1 - E0) / self.ħ, (E2 - E0) / self.ħ]

        delta = []
        J_sw  = []
        for n in (0, 1):
            D = [wc_vals[k] - omega[n] for k in (0, 1)]
            δ = sum(g_vals[k]**2 / D[k] for k in (0, 1))
            j = 0.5 * sum(g_vals[k] for k in (0, 1)) * sum(1 / D[k] for k in (0, 1))
            J_eff = j + self.g_cav
            delta.append(δ)
            J_sw.append(J_eff)

        H_mat = np.zeros((5, 5), dtype=complex)
        H_mat[0, 0] = 0.0
        H_mat[1, 1] = E1 + self.ħ * delta[0]
        H_mat[2, 2] = E1 + self.ħ * delta[0]
        H_mat[3, 3] = E2 + self.ħ * delta[1]
        H_mat[4, 4] = E2 + self.ħ * delta[1]
        H_mat[1, 2] = self.ħ * J_sw[0]; H_mat[2, 1] = self.ħ * J_sw[0]
        H_mat[3, 4] = self.ħ * J_sw[1]; H_mat[4, 3] = self.ħ * J_sw[1]

        return qt.Qobj(H_mat, dims=[[5], [5]])


# === 第二部分：LindbladSolver ===
class LindbladSolver:
    def __init__(self,
                 H: np.ndarray,
                 L: List[np.ndarray],
                 rho0: np.ndarray,
                 num_trajectories: int = 200,
                 dt: float = 0.01,
                 t_max: float = 5.0):
        """
        初始化 Lindblad 主方程 Monte Carlo 跳跃求解器
        :param H: Hamiltonian (NumPy 矩阵)
        :param L: 跳跃算符列表 (NumPy 矩阵)
        :param rho0: 初始密度矩阵
        """
        self.H = H
        self.L = L
        self.rho0 = rho0
        self.num_trajectories = num_trajectories
        self.dt = dt
        self.t_max = t_max
        self.t_points = np.arange(0, t_max, dt)
        self.dim = H.shape[0]
        self.state_densities = np.zeros((len(self.t_points), self.dim), dtype=float)

    def lindblad_rhs(self, rho: np.ndarray) -> np.ndarray:
        """
        计算 RHS: commutator + Lindblad 项
        """
        comm = -1j * (self.H @ rho - rho @ self.H)
        lind = np.zeros_like(rho, dtype=complex)
        for Lk in self.L:
            lind += (
                Lk @ rho @ Lk.conj().T
                - 0.5 * (Lk.conj().T @ Lk @ rho + rho @ Lk.conj().T @ Lk)
            )
        return comm + lind

    def monte_carlo_wavefunction(self, rho: np.ndarray) -> np.ndarray:
        """
        Monte Carlo 波函数跳跃：根据概率执行跳跃操作
        """
        jump_probs = np.array([
            np.real(np.trace(Lk @ rho @ Lk.conj().T)) for Lk in self.L
        ])
        total_jump_prob = jump_probs.sum()
        if total_jump_prob > 0 and np.random.rand() < total_jump_prob * self.dt:
            idx = np.random.choice(len(self.L), p=jump_probs / total_jump_prob)
            rho = self.L[idx] @ rho @ self.L[idx].conj().T
            rho = rho / np.trace(rho)
        return rho

    def solve(self):
        """
        求解并记录各态占据概率
        """
        for ti, _ in enumerate(self.t_points):
            rho_accum = np.zeros_like(self.rho0, dtype=complex)
            for _ in range(self.num_trajectories):
                rho = self.rho0.copy()
                for _ in self.t_points:
                    rho = rho + self.lindblad_rhs(rho) * self.dt
                    rho = self.monte_carlo_wavefunction(rho)
                rho_accum += rho
            rho_avg = rho_accum / self.num_trajectories
            self.state_densities[ti] = np.real(np.diag(rho_avg))


# === 主程序完整补充 ===
if __name__ == "__main__":
    cav_cfgs = [
        {"wavelength": 780e-9, "cavity_length": 20e-6,
         "cavity_waist": 2e-6, "finesse": 1000, "gamma_sp": 2 * np.pi * 6e6},
        {"wavelength": 852e-9, "cavity_length": 25e-6,
         "cavity_waist": 3e-6, "finesse": 1500, "gamma_sp": 2 * np.pi * 5e6}
    ]
    atom_levels = [0.0, 2 * np.pi * 6e9, 2 * np.pi * 12e9]

    #设置腔-腔耦合强度
    g_cav = 1e5  
    system = TwoCavityThreeLevelAtom(cav_cfgs, atom_levels, cavity_coupling=g_cav)
    H_sw_qobj = system.generate_sw_effective_hamiltonian()
    H_sw = H_sw_qobj.full()

    #设置跳跃算符
    gamma1, gamma2 = 1e6, 0.5e6
    L_list = []
    L_list.append((np.sqrt(gamma1) * qt.basis(5, 0) * qt.basis(5, 1).dag()).full())
    L_list.append((np.sqrt(gamma1) * qt.basis(5, 0) * qt.basis(5, 2).dag()).full())
    L_list.append((np.sqrt(gamma2) * qt.basis(5, 0) * qt.basis(5, 3).dag()).full())
    L_list.append((np.sqrt(gamma2) * qt.basis(5, 0) * qt.basis(5, 4).dag()).full())

    rho0_qobj = qt.basis(5, 1) * qt.basis(5, 1).dag()
    rho0 = rho0_qobj.full()

    #创建LindbladSolver实例并求解
    solver = LindbladSolver(H=H_sw,
                            L=L_list,
                            rho0=rho0,
                            num_trajectories=200,
                            dt=1e-9,
                            t_max=1e-6)
    solver.solve()

    plt.figure(figsize=(10, 6))
    for idx in range(solver.dim):
        plt.plot(solver.t_points * 1e6, solver.state_densities[:, idx],
                 label=f"State {idx}")

    plt.xlabel("Time (µs)")
    plt.ylabel("Occupation Probability")
    plt.title("Occupation Probability Evolution of Atomic States")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
