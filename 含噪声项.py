import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Callable
import scipy.sparse as sp

# === TwoCavityThreeLevelAtom ===
class TwoCavityThreeLevelAtom:
    """
    双腔耦合三能级原子-QED模型（含多跃迁、双腔耦合、激光驱动、失谐、噪声及耗散）
      - 每个光学腔光子数 0/1 (fock_dim=2)
      - 原子 3 能级
      - 总 Hilbert 空间: 2_cavA × 2_cavB × 3_atom = 12 维
      - 参数: λ, L, w0, F, γ_sp → FSR, Δν, Vmod, g_max, ω_c, κ, χ
      - 考虑跃迁: |0⟩↔|1⟩, |1⟩↔|2⟩, |0⟩↔|2⟩
      - 含腔-腔耦合 g_cav
      - 激光驱动 Rabi 频率 Ω, 驱动频率 ω_L (时依赖)
      - 考虑腔失谐 Δcav, 原子失谐 Δ10, Δ21
      - 添加噪声模型:
          • 白噪声 (Ornstein–Uhlenbeck 过程)
          • 耦合误差 (随机 Telegraph 噪声)
          • 串扰误差 (1/f 频谱噪声)
          • 退相干误差 (多通道 Lindblad)
      - 含腔耗散 κ_A, κ_B 与原子自发辐射 γ_sp
    """
    c = 2.99792458e8      # 光速 (m/s)
    ħ = 1.0545718e-34     # 约化普朗克常数

    def __init__(
        self,
        cavity_configs: List[Dict[str, float]],
        atom_levels: List[float],
        cavity_coupling: float = 0.0,
        noise_params: Dict[str, float] = None
    ):
        assert len(cavity_configs) == 2, "需要两个腔的配置"
        assert len(atom_levels) == 3, "需要三能级原子能级"
        self.atom_levels = atom_levels
        self.g_cav = cavity_coupling
        gamma_sp_list = [cfg.get("gamma_sp", 0.0) for cfg in cavity_configs]
        self.gamma_sp = np.mean(gamma_sp_list)
        
        # 纯退相干率
        self.gamma_phi_10 = 1.0
        self.gamma_phi_21 = 1.0

        # 模式维度与算符
        self.d_cav = 2
        self.d_atom = 3
        self.a = qt.destroy(self.d_cav)
        self.I_cav = qt.qeye(self.d_cav)
        self.I_atom = qt.qeye(self.d_atom)

        # 原子投影与跃迁算符
        self.proj = [qt.basis(self.d_atom, i) * qt.basis(self.d_atom, i).dag()
                     for i in range(self.d_atom)]
        self.sigma10 = qt.basis(self.d_atom, 0) * qt.basis(self.d_atom, 1).dag()
        self.sigma21 = qt.basis(self.d_atom, 2) * qt.basis(self.d_atom, 1).dag()
        self.sigma20 = qt.basis(self.d_atom, 2) * qt.basis(self.d_atom, 0).dag()

        # 腔参数列表
        self.params = []
        for cfg in cavity_configs:
            lam = cfg["wavelength"]
            L = cfg["cavity_length"]
            w0 = cfg["cavity_waist"]
            F = cfg["finesse"]
            γ_sp = cfg.get("gamma_sp", 0.0)
            fsr = self.c / (2 * L)
            Δν = fsr / F
            Vmod = np.pi * w0**2 * L / 4
            g_max = np.sqrt(3 * lam**2 * γ_sp / (8 * np.pi * Vmod))
            χ = cfg.get("chi", 0.0)  # Kerr 非线性系数
            ω_c = 2 * np.pi * (self.c / lam)
            κ = np.pi * Δν
            self.params.append({
                "wavelength": lam,
                "cavity_length": L,
                "cavity_waist": w0,
                "finesse": F,
                "FSR": fsr,
                "linewidth": Δν,
                "mode_volume": Vmod,
                "g_max": g_max,
                "chi": χ,
                "ω_c": ω_c,
                "κ": κ
            })
        
        # 噪声参数
        noise_params = noise_params or {}
        
        # OU 白噪声参数
        self.ou_tau = noise_params.get('ou_tau', 1e-6)
        self.ou_sigma = noise_params.get('white_noise_amp', 0.0)
        self.noise_dt = noise_params.get('noise_dt', 1e-9)
        self.ou_x1 = 0.0
        self.ou_x2 = 0.0
        self.coupling_err = 1.0
        
        # Telegraph 噪声参数
        self.telegraph_rate = noise_params.get('telegraph_rate', 1e3)
        self.telegraph_state = 1
        
        # 串扰 1/f 噪声参数
        self.crosstalk_err = noise_params.get('crosstalk_error', 0.0)
        self.crosstalk_phis = np.random.uniform(0, 2*np.pi, 4)
        self.crosstalk_omegas = np.logspace(3, 6, 4)
        
        # 退相干误差
        self.decohere_err_10 = noise_params.get('decoherence_10', 0.0)
        self.decohere_err_21 = noise_params.get('decoherence_21', 0.0)

        # 模式维度与算符
        self.d_cav, self.d_atom = 2, 3
        self.a = qt.destroy(self.d_cav)
        self.I_cav, self.I_atom = qt.qeye(self.d_cav), qt.qeye(self.d_atom)
        self.proj = [qt.basis(self.d_atom,i)*qt.basis(self.d_atom,i).dag() for i in range(self.d_atom)]
        self.sigma10 = qt.basis(self.d_atom,0)*qt.basis(self.d_atom,1).dag()
        self.sigma21 = qt.basis(self.d_atom,2)*qt.basis(self.d_atom,1).dag()
        self.sigma20 = qt.basis(self.d_atom,2)*qt.basis(self.d_atom,0).dag()

        # 腔参数
        self.params = []
        for cfg in cavity_configs:
            lam, L, w0 = cfg['wavelength'], cfg['cavity_length'], cfg['cavity_waist']
            F, γ_sp = cfg['finesse'], cfg.get('gamma_sp',0.0)
            fsr, Δν = self.c/(2*L), fsr/F
            Vmod = np.pi*w0**2*L/4
            g_max = np.sqrt(3*lam**2*γ_sp/(8*np.pi*Vmod))
            χ = cfg.get('chi',0.0)
            ω_c, κ = 2*np.pi*(self.c/lam), np.pi*Δν
            self.params.append(dict(
                wavelength=lam, cavity_length=L, cavity_waist=w0,
                finesse=F, FSR=fsr, linewidth=Δν, mode_volume=Vmod,
                g_max=g_max, chi=χ, ω_c=ω_c, κ=κ
            ))

    @staticmethod
    def get_cavity_operators() -> Tuple[qt.Qobj, qt.Qobj]:
        """
        返回截断到fock_dim=2子空间的腔湮灭和产生算符(a, a†)。
        """
        a_mat = np.array([[0,1],[0,0]])
        return qt.Qobj(a_mat), qt.Qobj(a_mat.T)

    @staticmethod
    def compute_photon_overlap(waveA: np.ndarray, waveB: np.ndarray, t: np.ndarray) -> float:
        """
        计算波包重叠度 γ = ∫ ψ_A*(t) ψ_B(t) dt，用于后选保真度统计。
        """
        dt = t[1]-t[0]
        return float(np.sum(np.conj(waveA)*waveB)*dt)

    @staticmethod
    def generate_single_photon_wavefunction(
        t_range: np.ndarray,
        shape: str = 'exponential',
        t0: float = 0.0,
        tau: float = None,
        sigma: float = None
    ) -> np.ndarray:
        """
        生成单光子时域波包 ψ(t)：
        - 指数型: \(ψ_{exp}(t)=Θ(t−t0)\sqrt{2/τ} e^{-(t−t0)/τ}\) (τ=波包长度≈1/κ) citeturn5file17
        - 高斯型: \(ψ_{gauss}(t)=1/(2πσ^2)^{1/4}\exp[-(t−t0)^2/(4σ^2)]\) (FWHM≈2√{2ln2}σ)
        参数:
          t_range: 时间网格
          shape: 'exponential' 或 'gaussian'
          t0: 波包起始时间
          tau: 指数衰减常数 (必需 for 'exponential')
          sigma: 高斯宽度 (必需 for 'gaussian')
        返回: 归一化波包数组
        """
        psi = np.zeros_like(t_range, dtype=complex)
        if shape == 'exponential':
            assert tau is not None, '需要指定 tau'
            norm = np.sqrt(2.0/tau)
            psi = np.where(
                t_range>=t0,
                norm * np.exp(-(t_range - t0)/tau),
                0.0
            )
        elif shape == 'gaussian':
            assert sigma is not None, '需要指定 sigma'
            norm = 1.0/(2*np.pi*sigma**2)**0.25
            psi = norm * np.exp(-(t_range - t0)**2/(4*sigma**2))
        else:
            raise ValueError('不支持的 shape 类型')
        # 归一化
        dt = t_range[1] - t_range[0]
        psi /= np.sqrt(np.sum(np.abs(psi)**2)*dt)
        return psi

    def ornstein_uhlenbeck_noise(self):
        """更新 OU 噪声过程"""
        θ, σ, dt = 1/self.ou_tau, self.ou_sigma, self.noise_dt
        e = np.exp(-θ*dt)
        ξ = np.random.normal(0,1)
        self.ou_x1 = e*self.ou_x1 + σ*np.sqrt(1-e**2)*ξ
        self.ou_x2 = e*self.ou_x2 + σ*np.sqrt(1-e**2)*ξ

    def white_noise_hamiltonian(self, t: float) -> qt.Qobj:
        """带相关性的 OU 噪声"""
        self.ornstein_uhlenbeck_noise()
        a1 = qt.tensor(self.a, self.I_cav, self.I_atom)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        return (self.ou_x1*(a1+a1.dag()) + self.ou_x2*(a2+a2.dag()))

    def coupling_error_hamiltonian(self) -> qt.Qobj:
        """Telegraph 随机翻转腔-腔耦合"""
        p = 1 - np.exp(-self.telegraph_rate*self.noise_dt)
        if np.random.rand() < p:
            self.telegraph_state *= -1
        a1 = qt.tensor(self.a, self.I_cav, self.I_atom)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        return (self.telegraph_state*self.coupling_err)*(a1.dag()*a2 + a1*a2.dag())

    def crosstalk_error_hamiltonian(self, t: float) -> qt.Qobj:
        """1/f 噪声模拟串扰"""
        total = 0
        for amp, ω, φ in zip([self.crosstalk_err]*4, self.crosstalk_omegas, self.crosstalk_phis):
            total += amp*(np.sin(ω*t+φ)/ω)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        s10 = qt.tensor(self.I_cav, self.I_cav, self.sigma10)
        return total*(a2.dag()*s10 + a2*s10.dag())

    def decoherence_error_operators(self) -> List[qt.Qobj]:
        """多通道非 Markovian 退相干"""
        ops = []
        # 10 跃迁退相干
        if self.decohere_err_10 > 0:
            rate = self.decohere_err_10*(1+0.1*np.random.randn())
            ops.append(np.sqrt(rate)*self.sigma10.dag()*self.sigma10)
        # 21 跃迁退相干
        if self.decohere_err_21 > 0:
            rate = self.decohere_err_21*(1+0.1*np.random.randn())
            ops.append(np.sqrt(rate)*self.sigma21.dag()*self.sigma21)
        return ops

    def laser_drive(self, Omega: float, omega_L: float, t: float) -> qt.Qobj:
        """
        时依赖激光驱动: H_drive(t) = ħ Ω[a e^{iω_L t}+a† e^{-iω_L t}]⊗I⊗I
        """
        return self.ħ*Omega*(
            qt.tensor(self.a, self.I_cav, self.I_atom)*np.exp(1j*omega_L*t) +
            qt.tensor(self.a.dag(), self.I_cav, self.I_atom)*np.exp(-1j*omega_L*t)
        )

    def total_hamiltonian(
        self,
        Omega: float,
        omega_L: float,
        t: float,
        g_max_10: float,
        g_max_21: float = 0.0,
        g_max_20: float = 0.0
    ) -> qt.Qobj:
        """
        返回时刻 t 的全 Hamiltonian (12×12)
        H_tot = ħ(H_cav1+H_cav2+H_nl+H_atom+H_int+H_cav_cav) + H_drive(t) + H_decay
        """
        E0, E1, E2 = self.atom_levels
        H_atom = qt.tensor(
            self.I_cav, self.I_cav,
            E0*self.proj[0] + E1*self.proj[1] + E2*self.proj[2]
        )
        a1 = qt.tensor(self.a, self.I_cav, self.I_atom)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        H_cav1 = self.params[0]["ω_c"] * a1.dag()*a1
        H_cav2 = self.params[1]["ω_c"] * a2.dag()*a2
        # Kerr 非线性项
        chi1 = self.params[0]["chi"]
        chi2 = self.params[1]["chi"]
        H_nl = self.ħ * (
            chi1 * a1.dag()*a1.dag()*a1*a1 +
            chi2 * a2.dag()*a2.dag()*a2*a2
        )
        s10 = qt.tensor(self.I_cav, self.I_cav, self.sigma10)
        s21 = qt.tensor(self.I_cav, self.I_cav, self.sigma21)
        s20 = qt.tensor(self.I_cav, self.I_cav, self.sigma20)
        H_int = self.ħ*(
            g_max_10*(a1.dag()*s10 + a1*s10.dag()) +
            g_max_21*(a1.dag()*s21 + a1*s21.dag()) +
            g_max_20*(a1.dag()*s20 + a1*s20.dag())
        )
        H_cav_cav = self.ħ*self.g_cav*(a1.dag()*a2 + a1*a2.dag())
        H_drive = self.laser_drive(Omega, omega_L, t)
        H_decay = -1j*0.5*(
            self.params[0]["κ"]*a1.dag()*a1 +
            self.params[1]["κ"]*a2.dag()*a2
        )

        H_0 = self.ħ*(H_cav1+H_cav2) + H_nl + H_atom + H_int + H_cav_cav + H_drive + H_decay
        H_deph = -1j*0.5*(
            self.params[0]["κ"]*a1.dag()*a1 + self.params[1]["κ"]*a2.dag()*a2
        )
        H_noise = self.white_noise_hamiltonian(t) + self.coupling_error_hamiltonian() + self.crosstalk_error_hamiltonian(t)
        return H0 + H_deph + H_noise

    def generate_effective_hamiltonian(
        self,
        detune_cav1: float = 0.0,
        detune_cav2: float = 0.0,
        detune_10: float = 0.0,
        detune_21: float = 0.0,
        Omega_eff: float = 0.0,
        g10: float = None,
        g21: float = None,
        g20: float = None
    ) -> qt.Qobj:
        """
        构造静态有效 H_eff (12×12) 包括失谐和有效驱动：
        H_eff = ħ(H_cav1+H_cav2 + H_atom_rot + H_int + H_cav_cav + H_drive_eff)
                - iħ/2(κ1 a1†a1 + κ2 a2†a2 + γ_sp(σ11+σ22))
        参数:
          detune_cav1/2: 腔失谐 Δcav = ω_c-ω_L
          detune_10/21: 原子跃迁失谐
          Omega_eff: 驱动有效 Rabi 频率
          g10/g21/g20: 各跃迁耦合 (默认取 params 中 g_max)
        """
        E0, E1, E2 = self.atom_levels
        #旋转系数
        H_atom_rot = qt.tensor(
            self.I_cav, self.I_cav,
            detune_10*self.proj[1] + (detune_10+detune_21)*self.proj[2]
        )
        
        a1 = qt.tensor(self.a, self.I_cav, self.I_atom)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        H_cav1 = (self.params[0]["ω_c"]-detune_cav1)*a1.dag()*a1
        H_cav2 = (self.params[1]["ω_c"]-detune_cav2)*a2.dag()*a2
        
        # Kerr 非线性项
        chi1 = self.params[0]["chi"]
        chi2 = self.params[1]["chi"]
        H_nl = self.ħ * (
            chi1 * a1.dag()*a1.dag()*a1*a1 +
            chi2 * a2.dag()*a2.dag()*a2*a2
        )

        #跃迁算符
        s10 = qt.tensor(self.I_cav, self.I_cav, self.sigma10)
        s21 = qt.tensor(self.I_cav, self.I_cav, self.sigma21)
        s20 = qt.tensor(self.I_cav, self.I_cav, self.sigma20)

        #耦合因子
        g10 = g10 or self.params[0]["g_max"]
        g21 = g21 or self.params[0]["g_max"]
        g20 = g20 or self.params[0]["g_max"]
        H_int = self.ħ*(
            g10*(a1.dag()*s10 + a1*s10.dag()) +
            g21*(a1.dag()*s21 + a1*s21.dag()) +
            g20*(a1.dag()*s20 + a1*s20.dag())
        )

        #静态驱动 Stark 移位
        H_drive_eff = self.ħ*(Omega_eff/2)*(s10 + s10.dag())

        #腔-腔耦合
        H_cav_cav = self.ħ*self.g_cav*(a1.dag()*a2 + a1*a2.dag())
        H_decay = -1j*0.5*(
            self.params[0]["κ"]*a1.dag()*a1 +
            self.params[1]["κ"]*a2.dag()*a2
        )
        
        H_eff = self.ħ*(H_cav1+H_cav2) + H_nl + H_atom_rot + H_int + H_cav_cav + H_drive_eff + H_decay
        
        # 添加噪声
        H_eff += self.coupling_error_hamiltonian() + self.crosstalk_error_hamiltonian(0.0)
        return H_eff

    def get_collapse_operators(self) -> List[qt.Qobj]:
        """
        返回 Lindblad 耗散算符: 腔衰减及原子自发衰减与原子纯退相干
        """
        a1 = qt.tensor(self.a, self.I_cav, self.I_atom)
        a2 = qt.tensor(self.I_cav, self.a, self.I_atom)
        s10 = qt.tensor(self.I_cav, self.I_cav, self.sigma10)
        s21 = qt.tensor(self.I_cav, self.I_cav, self.sigma21)
        ops = [
            np.sqrt(self.params[0]["κ"]) * a1,
            np.sqrt(self.params[1]["κ"]) * a2,
            np.sqrt(self.gamma_sp) * s10,
            np.sqrt(self.gamma_sp) * s21
        ]
        
        # 纯退相干
        if self.gamma_phi_10 > 0.0:
            ops.append(np.sqrt(self.gamma_phi_10) * qt.tensor(self.I_cav, self.I_cav, self.proj[1]))
        if self.gamma_phi_21 > 0.0:
            ops.append(np.sqrt(self.gamma_phi_21) * qt.tensor(self.I_cav, self.I_cav, self.proj[2]))
            
        # 额外退相干误差
        ops.extend(self.decoherence_error_operators())
        return ops

    def generate_sw_effective_hamiltonian(self) -> qt.Qobj:
        """
            基于 Schrieffer–Wolff 变换，在零光子子空间构造 5×5 原子有效哈密顿量。

            基底顺序 (原子态耦合情况)：
              0: |g,0_A,0_B>
              1: |e1,0_A,0_B>_A?  单激发态关联跨腔耦合 A
              2: |e1,0_A,0_B>_B?  单激发态关联跨腔耦合 B
              3: |e2,0_A,0_B>_A?  双激发态关联跨腔耦合 A
              4: |e2,0_A,0_B>_B?  双激发态关联跨腔耦合 B

            有效哈密顿量元件:
              H_eff[0,0] = E0
              H_eff[1,1] = E1 + ħ δ_10
              H_eff[2,2] = E1 + ħ δ_10
              H_eff[3,3] = E2 + ħ δ_20
              H_eff[4,4] = E2 + ħ δ_20
              H_eff[1,2] = H_eff[2,1] = ħ J_eff_10
              H_eff[3,4] = H_eff[4,3] = ħ J_eff_20

            其中:
              δ_n0 = Σ_{k=A,B} g_k,n^2 / Δ_k,n，原子 n→0 与腔 k 失谐下的二阶 Stark 位移
              J_eff_n0 = ½ Σ_k g_k,n * Σ_k 1/Δ_k,n + g_cav，包含交换耦合和直接腔-腔耦合

            返回 qt.Qobj 对象，dims=[[5],[5]]。
        """
        E0, E1, E2 = self.atom_levels
        omega10 = (E1 - E0)/self.ħ
        omega20 = (E2 - E0)/self.ħ
        g_vals = [p["g_max"] for p in self.params]
        wc_vals = [p["ω_c"] for p in self.params]
        delta = []
        
        J_eff = []
        for omega_n in [omega10, omega20]:
            D = [wc_vals[k] - omega_n for k in range(2)]
            δ_n0 = sum(g_vals[k]**2 / D[k] for k in range(2))
            J_n = 0.5 * sum(g_vals) * sum(1.0/D[k] for k in range(2))
            J_eff.append(J_n + self.g_cav)
            delta.append(δ_n0)
            
        H_mat = np.zeros((5,5),dtype=complex)
        H_mat[0,0] = E0
        H_mat[1,1] = E1 + self.ħ*delta[0]
        H_mat[2,2] = E1 + self.ħ*delta[0]
        H_mat[3,3] = E2 + self.ħ*delta[1]
        H_mat[4,4] = E2 + self.ħ*delta[1]
        H_mat[1,2] = H_mat[2,1] = self.ħ*J_eff[0]
        H_mat[3,4] = H_mat[4,3] = self.ħ*J_eff[1]
        return qt.Qobj(H_mat, dims=[[5],[5]])

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
            
# ===主程序===
if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import qutip as qt

    cav_cfgs = [
        {"wavelength": 780e-9, "cavity_length": 20e-6,
         "cavity_waist": 2e-6, "finesse": 1000, "gamma_sp": 2*np.pi*6e6, "chi": 0.0},
        {"wavelength": 852e-9, "cavity_length": 25e-6,
         "cavity_waist": 3e-6, "finesse": 1500, "gamma_sp": 2*np.pi*5e6, "chi": 0.0}
    ]
    atom_levels = [0.0, 2*np.pi*6e9, 2*np.pi*12e9]

    #设置腔-腔耦合强度
    g_cav = 1e5  
    system = TwoCavityThreeLevelAtom(cav_cfgs, atom_levels, cavity_coupling=g_cav)

    # 构造静态有效 Hamiltonian (12×12)
    H_eff_qobj = system.generate_effective_hamiltonian()
    H_eff = H_eff_qobj.full()

    # 构造 Schrieffer–Wolff 有效 Hamiltonian (5×5)
    H_sw_qobj = system.generate_sw_effective_hamiltonian()
    H_sw = H_sw_qobj.full()

    # 设置跳跃算符
    gamma1, gamma2 = 1e6, 0.5e6
    L_list = [
        (np.sqrt(gamma1) * qt.basis(5, 0) * qt.basis(5, 1).dag()).full(),
        (np.sqrt(gamma1) * qt.basis(5, 0) * qt.basis(5, 2).dag()).full(),
        (np.sqrt(gamma2) * qt.basis(5, 0) * qt.basis(5, 3).dag()).full(),
        (np.sqrt(gamma2) * qt.basis(5, 0) * qt.basis(5, 4).dag()).full()
    ]

    #初始密度矩阵：5维系统的第二个能级（索引 1）
    rho0_qobj = qt.basis(5, 1) * qt.basis(5, 1).dag()
    rho0 = rho0_qobj.full()

    #创建 LindbladSolver实例并求解（以 SW 哈密顿量为输入）
    solver = LindbladSolver(
        H=H_sw,
        L=L_list,
        rho0=rho0,
        num_trajectories=200,
        dt=1e-9,
        t_max=1e-6
    )
    solver.solve()

    plt.figure(figsize=(10, 6))
    for idx in range(solver.dim):
        plt.plot(
            solver.t_points * 1e6,
            solver.state_densities[:, idx],
            label=f"Level {idx}"
        )
    plt.xlabel("Time (µs)")
    plt.ylabel("Occupation Probability")
    plt.title("Occupation Probability Evolution (SW 5-dim system)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
