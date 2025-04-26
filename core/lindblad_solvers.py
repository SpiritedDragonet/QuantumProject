import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

class LindbladSolver:
    def __init__(self, H, L, rho0, dt=0.01, t_max=10, sparse=False):
        """
        初始化Lindblad主方程求解器
        :param H: 哈密顿量矩阵（二维NumPy数组）
        :param L: 跳跃算符的列表（每个元素是二维NumPy数组）
        :param rho0: 初始密度矩阵（二维NumPy数组）
        :param dt: 时间步长
        :param t_max: 最大时间
        :param sparse: 是否使用稀疏矩阵优化（用于大规模系统）
        """
        self.H = H
        self.L = L
        self.rho0 = rho0
        self.dt = dt
        self.t_max = t_max
        self.sparse = sparse
        self.dim = H.shape[0]
        self.t_points = np.arange(0, t_max, dt)
        self.rho = rho0
        self.results = []

        if sparse:
            self.H_sparse = sp.csr_matrix(H)
            self.L_sparse = [sp.csr_matrix(L_k) for L_k in L]
        else:
            self.H_sparse = H
            self.L_sparse = L

    def lindblad_master_equation(self, rho):
        """
        计算Lindblad主方程右边项，即密度矩阵的时间导数
        :param rho: 当前密度矩阵
        :return: d(rho)/dt
        """
        commutator = -1j * (np.dot(self.H, rho) - np.dot(rho, self.H))
        lindblad_term = np.zeros_like(rho, dtype=complex)
        for L_k in self.L_sparse:
            term = np.dot(L_k, np.dot(rho, L_k.conj().T)) - 0.5 * np.dot(np.dot(L_k.conj().T, L_k), rho) - 0.5 * np.dot(rho, np.dot(L_k.conj().T, L_k))
            lindblad_term += term

        return commutator + lindblad_term

    def monte_carlo_wavefunction(self, rho):
        """
        使用Monte Carlo波函数方法（MCWF）来模拟量子跳跃过程
        :param rho: 当前密度矩阵
        :return: 更新后的密度矩阵
        该函数的跃迁概率暂时有问题，处于NAN，调试中
        """
        jump_probs = np.array([np.sum(np.abs(np.dot(L_k, rho))**2) for L_k in self.L_sparse])
        total_jump_prob = np.sum(jump_probs)
        
        if total_jump_prob == 0:
            return rho  
            
        jump_event = np.random.rand() < total_jump_prob * self.dt

        if jump_event:
            jump_probs_normalized = jump_probs / total_jump_prob

            if np.any(np.isnan(jump_probs_normalized)) or np.any(np.isinf(jump_probs_normalized)):
                print("Warning: Invalid probability detected. Skipping jump.")
                return rho  

            chosen_index = np.random.choice(len(self.L_sparse), p=jump_probs_normalized)
            L_k = self.L_sparse[chosen_index]
            jump_strength = np.random.normal(1.0, 0.1)  

            rho = np.dot(L_k, np.dot(rho, L_k.conj().T)) * jump_strength     
        return rho

    def solve(self):
        """
        求解Lindblad主方程，并记录每个时间点的密度矩阵
        """
        for t in self.t_points:
            self.rho = self.rho + self.lindblad_master_equation(self.rho) * self.dt
            self.rho = self.monte_carlo_wavefunction(self.rho)
            self.results.append(np.real(np.trace(self.rho))) 

    def plot_results(self):
        """
        可视化求解结果（期望值随时间变化）
        """
        plt.plot(self.t_points, self.results)
        plt.xlabel('Time')
        plt.ylabel('Trace of Density Matrix')
        plt.title('Time Evolution of Expectation Value')
        plt.show()
