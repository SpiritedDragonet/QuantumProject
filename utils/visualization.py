"""
visualization_ui.py - 可视化UI模块，用于显示模拟结果
"""

import numpy as np
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider, Button, CheckButtons, TextBox
from mpl_toolkits.mplot3d import Axes3D
import time
import os
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import math
import glob
import h5py
from core.quantum_state import Frame

# 确保中文显示正常
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

class QuantumVisualization:
    """量子系统可视化类"""
    
    def __init__(self, config=None):
        """初始化可视化器"""
        self.config = config or {}
        
        # 创建结果目录
        self.results_path = self.config.get('save_path', './results')
        os.makedirs(self.results_path, exist_ok=True)
        
        # 动画参数
        self.animation = None
        self.paused = False
        
        # 保存的数据
        self.times = []
        self.states = []
        self.observables = {}

class ClassicalVisualization:
    """经典系统可视化UI类"""
    
    def __init__(self, sim_data=None):
        """初始化可视化UI
        
        参数:
            sim_data: SimulationData对象，用于UI和模拟之间的通信
        """
        self.sim_data = sim_data
        
        # 创建结果目录
        self.results_path = './results'
        os.makedirs(self.results_path, exist_ok=True)
        
        # 图形对象
        self.fig = None
        self.fig_3d = None
        self.ax_3d = None
        self.ax_energy = None
        self.ax_temp = None
        self.ax_info = None
        
        # 绘图元素
        self.trajectory_line = None
        self.current_point = None
        self.energy_line = None
        self.temp_line = None
        self.info_text = None
        
        # 控制元素
        self.pause_btn = None
        self.reset_btn = None
        self.cooling_btn = None
        self.history_btn = None
        self.window_size_var = None
        self.depth_var = None
        self.waist_var = None
        self.temp_var = None
        
        # 显示控制
        self.show_full_history = False  # 默认使用滚动窗口
        self.window_size_ms = 1000      # 默认窗口大小为1秒
        
        # 更新控制
        self.last_update_time = time.time()
        self.update_failures = 0
        
        # 物理常数
        self.kb = 1.380649e-23  # 玻尔兹曼常数 (J/K)
        
        # 回调函数
        self.callbacks = {}
    
    def setup_embedded_ui(self, fig, parent_window, callbacks=None):
        """在已有的Figure对象上设置UI
        
        参数:
            fig: matplotlib Figure对象
            parent_window: Tkinter窗口对象
            callbacks: 回调函数字典，包含UI事件处理函数
        """
        try:
            print("开始配置UI - 设置图形")
            # 存储图形和窗口对象
            self.fig = fig
            self.parent_window = parent_window
            self.callbacks = callbacks or {}
            
            # 清除任何现有的子图
            fig.clear()
            
            # 设置子图布局 - 2x2布局
            grid = fig.add_gridspec(2, 2, height_ratios=[2, 1], width_ratios=[2, 1])
            
            # 创建3D图形
            print("创建3D图形")
            from mpl_toolkits.mplot3d import Axes3D
            
            # 创建一个新的Figure用于3D图
            self.fig_3d = plt.figure(figsize=(5, 4))
            self.ax_3d = self.fig_3d.add_subplot(111, projection='3d')
            self.ax_3d.set_title('原子3D轨迹')
            self.ax_3d.set_xlabel('X (μm)')
            self.ax_3d.set_ylabel('Y (μm)')
            self.ax_3d.set_zlabel('Z (μm)')
            
            # 设置显示范围
            max_range = 5  # 微米
            self.ax_3d.set_xlim(-max_range, max_range)
            self.ax_3d.set_ylim(-max_range, max_range)
            self.ax_3d.set_zlim(-max_range, max_range)
            
            # 添加轨迹和当前点
            self.trajectory_line, = self.ax_3d.plot([], [], [], 'b-', lw=1, alpha=0.7)
            self.current_point, = self.ax_3d.plot([], [], [], 'ro', markersize=8)
            
            # 创建一个新的Frame用于3D视图
            self.frame_3d = ttk.LabelFrame(parent_window, text="3D视图")
            self.frame_3d.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=5)
            
            # 将3D图形嵌入到新Frame中
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
            self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, master=self.frame_3d)
            self.canvas_3d.draw()
            self.canvas_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            
            # 能量图
            print("创建能量图")
            self.ax_energy = fig.add_subplot(grid[0, 0])
            self.ax_energy.set_title('能量随时间变化')
            self.ax_energy.set_xlabel('时间 (ms)')
            self.ax_energy.set_ylabel('能量 (μK)')
            self.ax_energy.grid(True)
            self.energy_line, = self.ax_energy.plot([], [], 'g-')
            
            # 温度图
            print("创建温度图")
            self.ax_temp = fig.add_subplot(grid[0, 1])
            self.ax_temp.set_title('温度随时间变化')
            self.ax_temp.set_xlabel('时间 (ms)')
            self.ax_temp.set_ylabel('温度 (μK)')
            self.ax_temp.grid(True)
            self.temp_line, = self.ax_temp.plot([], [], 'r-')
            
            # 信息文本
            print("创建信息面板")
            self.ax_info = fig.add_subplot(grid[1, :])
            self.ax_info.set_axis_off()
            self.info_text = self.ax_info.text(0.05, 0.95, "正在初始化...", 
                                        transform=self.ax_info.transAxes, 
                                        verticalalignment='top', fontsize=10)
            
            # 获取初始参数
            if self.sim_data:
                params = self.sim_data.get_parameters()
                # 绘制光学陷阱
                self.plot_optical_trap(params['beam_waist'], params['rayleigh_length'])
                
                # 设置初始信息
                current_state = self.sim_data.get_current_state()
                self.update_info_text(
                    current_state['position'], 
                    current_state['velocity'], 
                    current_state['energy'] / self.kb * 1e6,  # 转换为μK
                    current_state['temperature'] * 1e6,       # 转换为μK
                    params['trap_depth'] * 1000,              # 转换为mK
                    params['beam_waist'], 
                    params['trap_freq_kHz'], 
                    params['enable_cooling']
                )
            
            # 创建控制UI组件
            print("创建控制组件")
            self.create_control_widgets(parent_window)
            
            # 调整布局
            fig.tight_layout()
            
            # 强制绘制，确保UI显示
            fig.canvas.draw()
            self.canvas_3d.draw()
            
            print("UI配置完成")
            return fig
            
        except Exception as e:
            print(f"设置UI时出错: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def create_control_widgets(self, parent_window):
        """创建控制UI组件
        
        参数:
            parent_window: Tkinter窗口对象
        """
        # 创建控制面板框架
        control_frame = ttk.Frame(parent_window)
        control_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        
        # 按钮组：暂停/继续、重置、启用/禁用冷却
        btn_frame = ttk.Frame(control_frame)
        btn_frame.pack(side=tk.LEFT, padx=10)
        
        # 获取当前参数
        params = self.sim_data.get_parameters() if self.sim_data else {}
        
        # 暂停/继续按钮
        self.pause_btn = ttk.Button(btn_frame, text="暂停", 
                                   command=lambda: self.call_callback('toggle_pause'))
        self.pause_btn.pack(side=tk.LEFT, padx=5)
        
        # 重置按钮
        self.reset_btn = ttk.Button(btn_frame, text="重置", 
                                   command=lambda: self.call_callback('reset_simulation'))
        self.reset_btn.pack(side=tk.LEFT, padx=5)
        
        # 冷却按钮
        enable_cooling = params.get('enable_cooling', False)
        self.cooling_btn = ttk.Button(btn_frame, 
                                     text="禁用冷却" if enable_cooling else "启用冷却", 
                                     command=lambda: self.call_callback('toggle_cooling'))
        self.cooling_btn.pack(side=tk.LEFT, padx=5)
        
        # 显示模式切换按钮
        self.history_btn = ttk.Button(btn_frame, 
                                     text="显示全部历史" if not self.show_full_history else "使用滚动窗口", 
                                     command=self.toggle_history_mode)
        self.history_btn.pack(side=tk.LEFT, padx=5)
        
        # 窗口大小调整
        window_frame = ttk.Frame(btn_frame)
        window_frame.pack(side=tk.LEFT, padx=5)
        ttk.Label(window_frame, text="窗口大小(ms):").pack(side=tk.LEFT)
        self.window_size_var = tk.StringVar(value=str(self.window_size_ms))
        window_entry = ttk.Entry(window_frame, textvariable=self.window_size_var, width=5)
        window_entry.pack(side=tk.LEFT, padx=2)
        window_entry.bind("<Return>", self.update_window_size)
        
        # 参数调整组
        param_frame = ttk.Frame(control_frame)
        param_frame.pack(side=tk.LEFT, padx=10, fill=tk.X, expand=True)
        
        # 阱深调整
        depth_frame = ttk.Frame(param_frame)
        depth_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(depth_frame, text="阱深 (mK):").pack(side=tk.LEFT)
        self.depth_var = tk.StringVar(value=str(params.get('trap_depth', 1.0) * 1000))
        depth_entry = ttk.Entry(depth_frame, textvariable=self.depth_var, width=7)
        depth_entry.pack(side=tk.LEFT, padx=5)
        depth_entry.bind("<Return>", lambda e: self.call_callback('update_trap_depth', self.depth_var.get()))
        
        # 束腰调整
        waist_frame = ttk.Frame(param_frame)
        waist_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(waist_frame, text="束腰 (μm):").pack(side=tk.LEFT)
        self.waist_var = tk.StringVar(value=str(params.get('beam_waist', 1.0)))
        waist_entry = ttk.Entry(waist_frame, textvariable=self.waist_var, width=7)
        waist_entry.pack(side=tk.LEFT, padx=5)
        waist_entry.bind("<Return>", lambda e: self.call_callback('update_beam_waist', self.waist_var.get()))
        
        # 温度调整
        temp_frame = ttk.Frame(param_frame)
        temp_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(temp_frame, text="温度 (μK):").pack(side=tk.LEFT)
        self.temp_var = tk.StringVar(value=str(params.get('temperature', 100e-6) * 1e6))
        temp_entry = ttk.Entry(temp_frame, textvariable=self.temp_var, width=7)
        temp_entry.pack(side=tk.LEFT, padx=5)
        temp_entry.bind("<Return>", lambda e: self.call_callback('update_temperature', self.temp_var.get()))
    
    def call_callback(self, callback_name, *args):
        """调用回调函数
        
        参数:
            callback_name: 回调函数名称
            *args: 回调函数参数
        """
        if self.callbacks and callback_name in self.callbacks:
            try:
                return self.callbacks[callback_name](*args)
            except Exception as e:
                print(f"调用回调函数 {callback_name} 出错: {e}")
                import traceback
                traceback.print_exc()
        else:
            print(f"回调函数 {callback_name} 不存在")
    
    def toggle_history_mode(self):
        """切换历史数据显示模式"""
        self.show_full_history = not self.show_full_history
        self.history_btn.configure(text="显示全部历史" if not self.show_full_history else "使用滚动窗口")
        self.update_visualization()  # 立即更新图表
    
    def update_window_size(self, event=None):
        """更新滚动窗口大小"""
        try:
            new_size = float(self.window_size_var.get())
            if new_size > 0:
                self.window_size_ms = new_size
                print(f"窗口大小已更新为 {new_size} ms")
                self.update_visualization()  # 立即更新图表
            else:
                print("窗口大小必须大于0")
                self.window_size_var.set(str(self.window_size_ms))  # 恢复原值
        except ValueError:
            print(f"无法解析窗口大小: {self.window_size_var.get()}")
            self.window_size_var.set(str(self.window_size_ms))  # 恢复原值
    
    def update_visualization(self):
        """更新可视化"""
        if not self.sim_data or not hasattr(self, 'fig') or self.fig is None:
            return
        
        # 限制更新频率
        current_time = time.time()
        if current_time - self.last_update_time < 0.05:  # 最多每50ms更新一次
            return
            
        self.last_update_time = current_time
        
        try:
            # 获取当前状态和历史数据
            current_state = self.sim_data.get_current_state()
            if not current_state['data_ready']:
                return
                
            history_data = self.sim_data.get_history_data()
            params = self.sim_data.get_parameters()
            
            # 提取数据
            position = current_state['position']
            velocity = current_state['velocity']
            energy = current_state['energy']
            temperature = current_state['temperature']
            
            # 提取历史数据
            positions = history_data['position_history']
            times = history_data['time_history']
            energies = history_data['energy_history'] / self.kb * 1e6  # 转换为μK
            temperatures = history_data['temperature_history'] * 1e6  # 转换为μK
            
            # 更新3D轨迹图 - 只显示当前点
            if hasattr(self, 'current_point') and self.current_point is not None:
                self.current_point.set_data([position[0]], [position[1]])
                self.current_point.set_3d_properties([position[2]])
            
            # 确保有数据可以绘制
            if len(times) > 0:
                # 转换时间单位为毫秒
                times_ms = times * 1000  # 转换为ms
                
                # 找到有效数据点 - 排除NaN和无穷大
                valid_indices = np.where(
                    np.isfinite(times_ms) & 
                    np.isfinite(energies) & 
                    np.isfinite(temperatures)
                )[0]
                
                if len(valid_indices) > 0:
                    valid_times = times_ms[valid_indices]
                    valid_energies = energies[valid_indices]
                    valid_temperatures = temperatures[valid_indices]
                    
                    # 直接绘制有效数据
                    if hasattr(self, 'energy_line') and self.energy_line is not None:
                        self.energy_line.set_data(valid_times, valid_energies)
                    
                    # 调整y轴范围，确保数据可见
                    if hasattr(self, 'ax_energy') and self.ax_energy is not None:
                        e_min = np.min(valid_energies)
                        e_max = np.max(valid_energies)
                        if e_max > e_min:
                            # 增加上下边距
                            range_e = e_max - e_min
                            y_min = max(0, e_min - 0.05 * range_e)
                            y_max = e_max + 0.05 * range_e
                            self.ax_energy.set_ylim(y_min, y_max)
                        
                        # 确保x轴范围合理 - 只显示实际数据所在的范围
                        x_min = np.min(valid_times)
                        x_max = np.max(valid_times)
                        if x_max > x_min:
                            # 增加左右边距
                            range_x = x_max - x_min
                            padding = 0.05 * range_x
                            self.ax_energy.set_xlim(x_min - padding, x_max + padding)
                
                    # 温度图同理
                    if hasattr(self, 'temp_line') and self.temp_line is not None:
                        self.temp_line.set_data(valid_times, valid_temperatures)
                    
                    if hasattr(self, 'ax_temp') and self.ax_temp is not None:
                        t_min = np.min(valid_temperatures)
                        t_max = np.max(valid_temperatures)
                        if t_max > t_min:
                            # 增加上下边距
                            range_t = t_max - t_min
                            y_min = max(0, t_min - 0.05 * range_t)
                            y_max = t_max + 0.05 * range_t
                            self.ax_temp.set_ylim(y_min, y_max)
                            
                        # 确保x轴范围与能量图一致
                        if hasattr(self, 'ax_energy'):
                            xlim = self.ax_energy.get_xlim()
                            self.ax_temp.set_xlim(xlim)
                
                    # 打印调试信息，帮助理解图表状态
                    min_time = np.min(valid_times)
                    max_time = np.max(valid_times)
                    print(f"图表数据: 点数={len(valid_times)}, 时间范围=[{min_time:.2f}, {max_time:.2f}]ms")
                    print(f"能量范围=[{np.min(valid_energies):.2f}, {np.max(valid_energies):.2f}]μK")
                    print(f"温度范围=[{np.min(valid_temperatures):.2f}, {np.max(valid_temperatures):.2f}]μK")
            
            # 更新信息文本
            energy_uK = energy / self.kb * 1e6  # 转换为μK
            temp_uK = temperature * 1e6  # 转换为μK
            self.update_info_text(
                position, velocity, energy_uK, temp_uK,
                params['trap_depth'] * 1000, params['beam_waist'],
                params['trap_freq_kHz'], params['enable_cooling']
            )
            
            # 更新图形
            if hasattr(self, 'fig') and self.fig is not None:
                self.fig.canvas.draw_idle()
                
            if hasattr(self, 'canvas_3d') and self.canvas_3d is not None:
                self.canvas_3d.draw_idle()
            
        except Exception as e:
            print(f"更新可视化时出错: {e}")
            import traceback
            traceback.print_exc()
            
            # 记录失败次数
            self.update_failures += 1
            if self.update_failures > 10:
                print("可视化更新失败次数过多，停止更新")
                return False
        
        return True
   
    def update_info_text(self, position, velocity, energy_uK, temp_uK, trap_depth_mK, beam_waist, trap_freq_kHz, enable_cooling):
        """更新信息文本
        
        参数:
            position: 原子位置（微米）
            velocity: 原子速度（米/秒）
            energy_uK: 能量（微开）
            temp_uK: 温度（微开）
            trap_depth_mK: 阱深（毫开）
            beam_waist: 束腰半径（微米）
            trap_freq_kHz: 陷阱频率（千赫兹）
            enable_cooling: 是否启用冷却
        """
        if not hasattr(self, 'info_text') or self.info_text is None:
            return
            
        self.info_text.set_text(
            f"位置: ({position[0]:.2f}, {position[1]:.2f}, {position[2]:.2f}) μm\n"
            f"速度: ({velocity[0]:.2f}, {velocity[1]:.2f}, {velocity[2]:.2f}) m/s\n"
            f"能量: {energy_uK:.2f} μK\n"
            f"温度: {temp_uK:.2f} μK\n"
            f"阱深: {trap_depth_mK:.2f} mK\n"
            f"束腰: {beam_waist:.2f} μm\n"
            f"频率: {trap_freq_kHz:.2f} kHz\n"
            f"冷却: {'开启' if enable_cooling else '关闭'}"
        )
    
    def plot_optical_trap(self, beam_waist, rayleigh_length):
        """绘制光学陷阱形状"""
        if not hasattr(self, 'ax_3d') or self.ax_3d is None:
            print("3D轴未创建，跳过光学陷阱绘制")
            return
        
        try:
            print(f"绘制光学陷阱: 束腰={beam_waist}μm, 瑞利长度={rayleigh_length*1e6:.2f}μm")
            
            # 清除原有的光阱可视化（保留轨迹和当前点）
            artists_to_keep = []
            if hasattr(self, 'trajectory_line'):
                artists_to_keep.append(self.trajectory_line)
            if hasattr(self, 'current_point'):
                artists_to_keep.append(self.current_point)
                
            # 安全移除图形元素
            for collection in list(self.ax_3d.collections):
                if collection not in artists_to_keep:
                    try:
                        collection.remove()
                    except:
                        pass
            
            for line in list(self.ax_3d.lines):
                if line not in artists_to_keep:
                    try:
                        line.remove()
                    except:
                        pass
            
            # 束腰面
            u = np.linspace(0, 2*np.pi, 20)
            v = np.linspace(0, np.pi, 20)
            x = beam_waist * np.outer(np.cos(u), np.sin(v))
            y = beam_waist * np.outer(np.sin(u), np.sin(v))
            z = np.zeros_like(x)
            
            self.ax_3d.plot_surface(x, y, z, color='r', alpha=0.1)
            
            # 光束轮廓
            z_range = np.linspace(-3*beam_waist, 3*beam_waist, 30)
            for zi in z_range[::5]:  # 每隔5个点画一个圆
                r = beam_waist * np.sqrt(1 + (zi*1e-6/rayleigh_length)**2)
                theta = np.linspace(0, 2*np.pi, 20)
                x_circle = r * np.cos(theta)
                y_circle = r * np.sin(theta)
                z_circle = zi * np.ones_like(theta)
                self.ax_3d.plot(x_circle, y_circle, z_circle, 'r-', alpha=0.3)
            
            # 中轴线
            self.ax_3d.plot([0, 0], [0, 0], [-3*beam_waist, 3*beam_waist], 'r--', alpha=0.3)
            
            # 更新3D图形
            if hasattr(self, 'canvas_3d') and self.canvas_3d is not None:
                try:
                    self.canvas_3d.draw_idle()
                except Exception as e:
                    print(f"更新3D画布失败: {e}")
                    
        except Exception as e:
            print(f"绘制光学陷阱时出错: {e}")
            import traceback
            traceback.print_exc()
    
    def update_button_states(self, paused=None, enable_cooling=None):
        """更新按钮状态
        
        参数:
            paused: 是否暂停
            enable_cooling: 是否启用冷却
        """
        if paused is not None and hasattr(self, 'pause_btn') and self.pause_btn:
            self.pause_btn.configure(text="继续" if paused else "暂停")
        
        if enable_cooling is not None and hasattr(self, 'cooling_btn') and self.cooling_btn:
            self.cooling_btn.configure(text="禁用冷却" if enable_cooling else "启用冷却")

# 添加DLCZ协议结果可视化函数
def plot_results(times, results, system_state, config, save_path=None):
    """
    绘制模拟结果图表
    
    参数：
        times (list): 时间点列表
        results (dict): 模拟结果字典
        system_state (SystemState): 系统状态对象
        config (Config): 配置对象
        save_path (str, optional): 图像保存路径
    """
    # 创建图形和子图
    fig = plt.figure(figsize=(14, 10))
    gs = plt.GridSpec(3, 2, figure=fig)
    
    # 1. 原子激发图
    ax1 = fig.add_subplot(gs[0, 0])
    for i, excitations in enumerate(results['atom_excitations']):
        ax1.plot(times * 1e6, excitations, label=f'原子 {i+1}')
    ax1.set_xlabel('时间 (μs)')
    ax1.set_ylabel('激发态布居数')
    ax1.set_title('原子激发动态')
    ax1.legend()
    ax1.grid(True)
    
    # 2. 腔光子数图
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(times * 1e6, results['cavity_photons'], 'r-')
    ax2.set_xlabel('时间 (μs)')
    ax2.set_ylabel('光子数')
    ax2.set_title('腔光子数演化')
    ax2.grid(True)
    
    # 3. 信号光子和闲置光子图
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(times * 1e6, results['signal_photons'], 'g-', label='信号光子')
    ax3.plot(times * 1e6, results['idler_photons'], 'b-', label='闲置光子')
    ax3.set_xlabel('时间 (μs)')
    ax3.set_ylabel('光子数')
    ax3.set_title('DLCZ协议光子生成')
    ax3.legend()
    ax3.grid(True)
    
    # 4. 纠缠保真度图
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(times * 1e6, results['entanglement_fidelity'], 'k-')
    ax4.set_xlabel('时间 (μs)')
    ax4.set_ylabel('保真度')
    ax4.set_title('量子纠缠保真度')
    ax4.set_ylim(0, 1.05)
    ax4.grid(True)
    
    # 5. 系统信息
    ax_info = fig.add_subplot(gs[2, :])
    ax_info.axis('off')
    
    # 添加系统信息文本
    info_text = (
        f"系统信息:\n"
        f"  原子类型: {config.atom_type}\n"
        f"  原子数量: {system_state.num_atoms}\n"
        f"  量子态维度: {system_state.system_dim}\n"
        f"  模拟类型: {config.simulation_type}\n"
        f"  模拟时长: {max(times)*1e6:.2f} μs\n"
    )
    if 'entanglement_fidelity' in results and len(results['entanglement_fidelity']) > 0:
        info_text += f"  最终纠缠保真度: {results['entanglement_fidelity'][-1]:.4f}\n"
    
    ax_info.text(0.01, 0.5, info_text, fontsize=10, verticalalignment='center')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图像
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_dlcz_protocol_results(protocol, figsize=(12, 10), save_path=None):
    """
    绘制DLCZ协议执行结果图表
    
    参数:
        protocol: DLCZProtocol对象，包含协议执行结果
        figsize: 图像大小
        save_path: 保存路径，如果不为None则保存图像
            
    返回:
        plt.Figure: 图像对象
    """
    fig, axes = plt.subplots(3, 1, figsize=figsize)
    
    # 1. 绘制写入过程结果
    if protocol.stage_results['write']:
        write_result = protocol.stage_results['write']
        times = write_result.get('times', [])
        
        if times and 'p_ground' in write_result and 'p_excited' in write_result and 'p_storage' in write_result:
            axes[0].plot(np.array(times)*1e6, write_result['p_ground'], 'b-', label='|g⟩')
            axes[0].plot(np.array(times)*1e6, write_result['p_excited'], 'r-', label='|e⟩')
            axes[0].plot(np.array(times)*1e6, write_result['p_storage'], 'g-', label='|s⟩')
            axes[0].set_ylabel('布居数')
            axes[0].set_title('写入过程')
            axes[0].legend()
            axes[0].grid(True, alpha=0.3)
            
            # 绘制信号光子
            if 'signal_photons' in write_result:
                ax_twin = axes[0].twinx()
                ax_twin.plot(np.array(times)*1e6, write_result['signal_photons'], 'm--', label='信号光子')
                ax_twin.set_ylabel('光子数')
                ax_twin.legend(loc='upper right')
    
    # 2. 绘制存储过程结果
    if protocol.stage_results['storage']:
        storage_result = protocol.stage_results['storage']
        times = storage_result.get('times', [])
        
        if times and 'p_storage' in storage_result:
            axes[1].plot(np.array(times)*1e6, storage_result['p_storage'], 'g-', label='|s⟩')
            axes[1].set_ylabel('存储态布居数')
            axes[1].set_title('存储过程')
            axes[1].legend()
            axes[1].grid(True, alpha=0.3)
            
            # 绘制相干性
            if 'coherence' in storage_result:
                ax_twin = axes[1].twinx()
                ax_twin.plot(np.array(times)*1e6, storage_result['coherence'], 'c--', label='相干性')
                ax_twin.set_ylabel('相干性')
                ax_twin.legend(loc='upper right')
    
    # 3. 绘制读取过程结果
    if protocol.stage_results['read']:
        read_result = protocol.stage_results['read']
        times = read_result.get('times', [])
        
        if times and 'p_ground' in read_result and 'p_excited' in read_result and 'p_storage' in read_result:
            axes[2].plot(np.array(times)*1e6, read_result['p_ground'], 'b-', label='|g⟩')
            axes[2].plot(np.array(times)*1e6, read_result['p_excited'], 'r-', label='|e⟩')
            axes[2].plot(np.array(times)*1e6, read_result['p_storage'], 'g-', label='|s⟩')
            axes[2].set_xlabel('时间 (μs)')
            axes[2].set_ylabel('布居数')
            axes[2].set_title('读取过程')
            axes[2].legend()
            axes[2].grid(True, alpha=0.3)
            
            # 绘制闲置光子
            if 'idler_photons' in read_result:
                ax_twin = axes[2].twinx()
                ax_twin.plot(np.array(times)*1e6, read_result['idler_photons'], 'm--', label='闲置光子')
                ax_twin.set_ylabel('光子数')
                ax_twin.legend(loc='upper right')
    
    # 添加纠缠验证结果作为标题
    if protocol.stage_results['verification']:
        verification = protocol.stage_results['verification']
        fidelity = verification.get('entanglement_fidelity', 0)
        bell_param = verification.get('bell_parameter', 0)
        success_prob = verification.get('success_probability', 0)
        
        fig.suptitle(f'DLCZ协议结果\n保真度: {fidelity:.4f}, Bell参数: {bell_param:.4f}, 成功概率: {success_prob:.1e}',
                   fontsize=14, y=0.98)
    else:
        fig.suptitle('DLCZ协议结果', fontsize=14, y=0.98)
    
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存图像
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
    return fig

def generate_call_graph(config):
    """
    使用pyan3自动生成调用关系图
    """
    try:
        print("正在使用pyan3生成调用关系图...")
        
        # 只获取项目的主要Python文件，排除外部库
        python_files = [
            'main.py',
            'config.py'
        ]
        
        # 添加各模块目录下的文件
        for module_dir in ['system', 'physics', 'dynamics', 'photonics', 'protocols', 'managers', 'utils']:
            if os.path.exists(module_dir):
                module_files = glob.glob(f'{module_dir}/*.py')
                python_files.extend(module_files)
        
        if not python_files:
            print("未找到Python文件，无法生成调用图")
            return
            
        print(f"将分析以下文件: {', '.join(python_files)}")
        
        # 使用pyan3生成调用关系图
        try:
            os.makedirs('debug_output', exist_ok=True)
            pyan3_cmd = [
                'pyan3',
                *python_files,
                '--uses',
                '--no-defines',
                '--colored',
                '--grouped',
                '--dot',
                '--output', 'debug_output/callgraph.dot'
            ]
            print(f"执行命令: {' '.join(pyan3_cmd)}")
            import subprocess
            subprocess.run(pyan3_cmd, check=True)
            print("使用pyan3生成的调用图已保存到debug_output/callgraph.dot")
        except Exception as e:
            print(f"使用pyan3失败: {e}")
            print("无法生成调用关系图")
    except Exception as e:
        print(f"生成调用关系图时出错: {e}")
        print("无法生成调用关系图")

def plot_excitation_probability(time_points, probs, title="原子激发概率", save_path=None):
    """
    绘制原子激发概率随时间的变化
    
    参数:
        time_points: 时间点数组
        probs: 激发概率数组，形状为[n, 3]，包含时间和两个原子的激发概率
        title: 图表标题
        save_path: 保存路径，如果为None则显示图表而不保存
        
    返回:
        matplotlib.figure.Figure: 图表对象
    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 绘制激发概率
    ax.plot(probs[:, 0], probs[:, 1], 'r-', label='B原子激发概率')
    ax.plot(probs[:, 0], probs[:, 2], 'b-', label='C原子激发概率')
    
    # 设置标签和标题
    ax.set_xlabel('时间 (ps)')
    ax.set_ylabel('激发概率')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 保存图表
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"图表已保存至 {save_path}")
    
    return fig

def plot_wavepackets(t_array, wavepackets, labels=None, title="波包可视化", figsize=(10, 6), save_path=None):
    """
    可视化多个波包在时间域上的分布
    
    参数:
        t_array (numpy.ndarray): 时间点数组
        wavepackets (list of numpy.ndarray): 波包数组列表，每个元素是一个波包的振幅
        labels (list of str, optional): 每个波包的标签
        title (str, optional): 图表标题
        figsize (tuple, optional): 图表大小
        save_path (str, optional): 如果提供，图表将保存到此路径
    
    返回:
        matplotlib.figure.Figure: 生成的图表对象
    """
    plt.figure(figsize=figsize)
    
    if labels is None:
        labels = [f"波包 {i+1}" for i in range(len(wavepackets))]
    
    for i, (wave, label) in enumerate(zip(wavepackets, labels)):
        # 绘制波包振幅
        plt.plot(t_array, np.abs(wave)**2, label=f"{label} 强度")
        
        # 添加波包相位的虚线表示
        if np.iscomplexobj(wave):
            phase = np.angle(wave)
            # 归一化相位以便图表展示
            phase_normalized = phase / np.pi
            plt.plot(t_array, phase_normalized, '--', label=f"{label} 相位/π")
    
    plt.xlabel("时间 (ns)")
    plt.ylabel("归一化振幅 / 相位")
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    
    # 添加波包特性的说明框
    textstr = "\n".join([
        f"波包 {i+1} 特性:" +
        f"\n  峰值位置: {t_array[np.argmax(np.abs(wave)**2)]:.2f} ns" +
        f"\n  有效宽度: {np.sum(np.abs(wave)**2) / np.max(np.abs(wave)**2) * (t_array[1] - t_array[0]):.2f} ns"
        for i, wave in enumerate(wavepackets)
    ])
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.annotate(textstr, xy=(0.05, 0.95), xycoords='axes fraction',
                 verticalalignment='top', horizontalalignment='left',
                 bbox=props, fontsize=9)
    
    # 如果提供了保存路径，则保存图表
    if save_path:
        try:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"波包图表已保存至: {save_path}")
        except Exception as e:
            print(f"保存波包图表时出错: {e}")
    
    return plt.gcf()

def create_density_matrix_plot(rho, ax=None, title="密度矩阵可视化", cmap='viridis', 
                              show_values=False, show_phase=True, show_magnitude=True,
                              system_labels=None, figsize=(10, 8), subsystems=None,
                              max_dim_to_show=30):
    """
    创建密度矩阵的可视化图
    
    参数:
        rho (numpy.ndarray 或 qutip.Qobj): 密度矩阵
        ax (matplotlib.axes.Axes, optional): 绘图使用的轴对象，默认为None（创建新图）
        title (str): 图表标题
        cmap (str): 颜色映射名称
        show_values (bool): 是否显示数值
        show_phase (bool): 是否显示相位图
        show_magnitude (bool): 是否显示幅值图
        system_labels (list): 系统标签，用于标注基矢
        figsize (tuple): 图像大小
        subsystems (tuple): 子系统维度列表，例如 (5, 3) 表示5维原子态和3维光子态的直积系统，总维度为15
        max_dim_to_show (int): 显示的最大维度，防止图像过大
        
    返回:
        matplotlib.figure.Figure: 图形对象
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    
    # 如果是qutip对象，转换为numpy数组
    if hasattr(rho, 'full'):
        rho_data = rho.full()
    else:
        rho_data = np.array(rho)
    
    # 处理密度矩阵维度
    dim = rho_data.shape[0]
    print(f"处理密度矩阵，维度: {dim}×{dim}")
    
    # 设置默认的子系统维度（如果未指定）
    if subsystems is None:
        # 根据密度矩阵维度推断子系统结构
        if dim == 15:
            subsystems = (5, 3)  # 5维原子态 × 3维光子态
            print("推断子系统为: 5维原子态 × 3维光子态")
        elif dim == 25:
            subsystems = (5, 5)  # 两个5维原子态
            print("推断子系统为: 5维原子态 × 5维原子态")
        elif dim == 225:  # 15×15直积系统
            subsystems = (15, 15)  # 两个15维子系统
            print("推断子系统为: 15维系统 × 15维系统")
        else:
            # 其他情况，假设是单一系统
            subsystems = (dim,)
            print(f"无法匹配已知模式，假设为单一系统: {dim}维")
            
    print(f"子系统维度设置为: {subsystems}，总维度: {dim}×{dim}")
    
    # 创建图形对象
    if ax is None:
        fig, axs = plt.subplots(1, 2 if show_phase else 1, figsize=figsize)
        if show_phase:
            ax_mag, ax_phase = axs if hasattr(axs, '__len__') else (axs, None)
        else:
            ax_mag = axs
            ax_phase = None
    else:
        # 确保ax是列表或元组
        if isinstance(ax, (list, tuple)):
            if len(ax) >= 2 and show_phase:
                ax_mag, ax_phase = ax[0], ax[1]
            else:
                ax_mag, ax_phase = ax[0], None
        else:
            ax_mag = ax
            ax_phase = None
        fig = ax_mag.figure
    
    # 计算幅值和相位
    magnitude = np.abs(rho_data)
    phase = np.angle(rho_data)
    
    # 创建基矢标签
    if system_labels is None:
        if len(subsystems) == 2:
            # 为两个子系统创建直积基矢标签
            d1, d2 = subsystems
            
            # 创建原子态标签 (将g3改为r表示里德堡态)
            if d1 == 5:
                atom_labels = ['$|g_1\\rangle$', '$|g_2\\rangle$', '$|e\\rangle$', '$|r\\rangle$', '$|g_3\\rangle$']
            else:
                atom_labels = [f'$|{i}\\rangle$' for i in range(d1)]
                
            # 创建光子态标签
            if d2 == 3:
                photon_labels = ['$|0\\rangle$', '$|1\\rangle$', '$|2\\rangle$']
            else:
                photon_labels = [f'$|{i}\\rangle$' for i in range(d2)]
            
            # 创建直积基矢标签
            labels = []
            for i in range(min(d1, 10)):  # 限制最多显示10个
                for j in range(min(d2, 10)):  # 限制最多显示10个
                    # 使用更清晰的组合标签格式
                    atom = atom_labels[i].replace('$|', '').replace('\\rangle$', '')
                    photon = photon_labels[j].replace('$|', '').replace('\\rangle$', '')
                    labels.append(f'$|{atom},{photon}\\rangle$')
            
            # 如果标签过多，限制显示数量
            if len(labels) > max_dim_to_show:
                print(f"标签数量 {len(labels)} 超过显示限制 {max_dim_to_show}，将只显示部分标签")
                step = max(1, len(labels) // max_dim_to_show)
                labels = labels[::step]
        else:
            # 单个系统的标签
            labels = [f'$|{i}\\rangle$' for i in range(min(dim, max_dim_to_show))]
    else:
        # 使用用户提供的标签
        labels = system_labels
    
    # 裁剪密度矩阵（如果需要）
    if dim > max_dim_to_show:
        print(f"矩阵尺寸 {dim}×{dim} 超过显示限制 {max_dim_to_show}，将裁剪显示")
        rho_data = rho_data[:max_dim_to_show, :max_dim_to_show]
        magnitude = magnitude[:max_dim_to_show, :max_dim_to_show]
        phase = phase[:max_dim_to_show, :max_dim_to_show]
        labels = labels[:max_dim_to_show]
        dim = max_dim_to_show
    
    # 绘制幅值图
    if show_magnitude:
        im_mag = ax_mag.imshow(magnitude, cmap=cmap, interpolation='nearest')
        ax_mag.set_title(f'{title} - 幅值')
        fig.colorbar(im_mag, ax=ax_mag)
        
        # 设置刻度标签（只在维度小于一定值时显示全部标签）
        if dim <= max_dim_to_show:
            ax_mag.set_xticks(np.arange(dim))
            ax_mag.set_yticks(np.arange(dim))
            ax_mag.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax_mag.set_yticklabels(labels, fontsize=8)
        else:
            # 对于大维度矩阵，只显示部分标签
            step = max(1, dim // 10)
            tick_positions = np.arange(0, dim, step)
            ax_mag.set_xticks(tick_positions)
            ax_mag.set_yticks(tick_positions)
            
            # 选择对应位置的标签
            tick_labels = [labels[i] if i < len(labels) else f'$|{i}\\rangle$' for i in tick_positions]
            ax_mag.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8)
            ax_mag.set_yticklabels(tick_labels, fontsize=8)
        
        # 在矩阵上标注数值
        if show_values and dim <= 10:  # 只在矩阵较小时显示数值
            for i in range(dim):
                for j in range(dim):
                    if magnitude[i, j] > 0.01:  # 只显示幅值大于阈值的元素
                        ax_mag.text(j, i, f'{magnitude[i, j]:.2f}',
                                   ha='center', va='center', color='w' if magnitude[i, j] > 0.5 else 'k')
    
    # 绘制相位图
    if show_phase and ax_phase is not None:
        # 使用循环色图显示相位
        phase_cmap = 'hsv'  # 循环色图适合显示角度
        im_phase = ax_phase.imshow(phase, cmap=phase_cmap, interpolation='nearest', vmin=-np.pi, vmax=np.pi)
        ax_phase.set_title(f'{title} - 相位')
        fig.colorbar(im_phase, ax=ax_phase)
        
        # 设置刻度标签
        if dim <= max_dim_to_show:
            ax_phase.set_xticks(np.arange(dim))
            ax_phase.set_yticks(np.arange(dim))
            ax_phase.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
            ax_phase.set_yticklabels(labels, fontsize=8)
        else:
            # 对于大维度矩阵，只显示部分标签
            step = max(1, dim // 10)
            tick_positions = np.arange(0, dim, step)
            ax_phase.set_xticks(tick_positions)
            ax_phase.set_yticks(tick_positions)
            
            # 选择对应位置的标签
            tick_labels = [labels[i] if i < len(labels) else f'$|{i}\\rangle$' for i in tick_positions]
            ax_phase.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8)
            ax_phase.set_yticklabels(tick_labels, fontsize=8)
    
    plt.tight_layout()
    return fig

def create_density_matrix_animation(evolution_data, time_points=None, title="密度矩阵时间演化", 
                                   subsystems=(15, 15), save_path=None, fps=5, dpi=100):
    """
    创建密度矩阵随时间演化的动画
    
    参数:
        evolution_data: 列表或字典，包含不同时间点的密度矩阵
        time_points: 时间点数组，长度应与evolution_data相同
        title: 动画标题
        subsystems: 子系统维度，例如(15, 15)表示两个15维系统
        save_path: 动画保存路径，如果为None则不保存
        fps: 帧率
        dpi: 图像分辨率
        
    返回:
        matplotlib.animation.FuncAnimation: 动画对象
    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import numpy as np
    
    try:
        # 创建独立的图形窗口，避免嵌入到Tkinter中导致的包装冲突
        fig = plt.figure(figsize=(14, 6))
        
        # 创建两个子图
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        
        # 获取第一帧数据用于初始化
        if isinstance(evolution_data, dict):
            # 如果是字典，按时间点排序
            times = sorted(evolution_data.keys())
            rho_first = evolution_data[times[0]]
        else:
            # 如果是列表，直接使用第一个元素
            rho_first = evolution_data[0]
            if time_points is None:
                times = np.arange(len(evolution_data))
            else:
                times = time_points
        
        # 初始化幅值和相位图
        if hasattr(rho_first, 'full'):
            rho_data = rho_first.full()
        else:
            rho_data = np.array(rho_first)
        
        magnitude = np.abs(rho_data)
        phase = np.angle(rho_data)
        
        im_mag = ax1.imshow(magnitude, cmap='viridis', interpolation='nearest')
        im_phase = ax2.imshow(phase, cmap='hsv', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
        
        colorbar_mag = plt.colorbar(im_mag, ax=ax1)
        colorbar_phase = plt.colorbar(im_phase, ax=ax2)
        
        ax1.set_title(f'{title} - 幅值')
        ax2.set_title(f'{title} - 相位')
        
        # 创建时间标题
        time_text = fig.suptitle(f'时间: {times[0]:.3e}', fontsize=16)
        
        # 设置标签（针对较小维度的矩阵）
        dim = rho_data.shape[0]
        if dim <= 30 and subsystems is not None:
            d1, d2 = subsystems
            labels = []
            # 创建简化标签，选择部分量子态展示
            step1 = max(1, d1 // 5)
            step2 = max(1, d2 // 5)
            for i in range(0, min(d1, 10), step1):
                for j in range(0, min(d2, 10), step2):
                    labels.append(f'$|{i},{j}\\rangle$')
            
            # 截断标签列表以匹配矩阵维度
            labels = labels[:dim]
            
            if len(labels) <= 30:
                for ax in [ax1, ax2]:
                    ax.set_xticks(np.arange(len(labels)))
                    ax.set_yticks(np.arange(len(labels)))
                    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
                    ax.set_yticklabels(labels, fontsize=8)
        
        # 更新函数
        def update(frame):
            if isinstance(evolution_data, dict):
                rho = evolution_data[times[frame]]
                t = times[frame]
            else:
                rho = evolution_data[frame]
                t = times[frame]
            
            if hasattr(rho, 'full'):
                rho_data = rho.full()
            else:
                rho_data = np.array(rho)
            
            magnitude = np.abs(rho_data)
            phase = np.angle(rho_data)
            
            im_mag.set_array(magnitude)
            im_phase.set_array(phase)
            time_text.set_text(f'时间: {t:.3e}')
            
            return im_mag, im_phase, time_text
        
        # 创建动画
        ani = animation.FuncAnimation(fig, update, frames=len(times), 
                                     interval=1000/fps, blit=True)
        
        plt.tight_layout()
        fig.subplots_adjust(top=0.9)  # 为标题留出空间
        
        # 保存动画
        if save_path:
            try:
                ani.save(save_path, fps=fps, dpi=dpi, writer='pillow')
                print(f"动画已保存至 {save_path}")
            except Exception as e:
                print(f"保存动画时出错: {e}")
                
        # 显示图形（无阻塞）
        plt.show(block=False)
        
        return ani
    except Exception as e:
        import traceback
        print(f"创建密度矩阵动画时出错: {e}")
        traceback.print_exc()
        return None

def visualize_density_matrix_evolution(frames, system_name=None, title="量子态演化", 
                                       subsystems=(15, 15), figsize=(12, 10), 
                                       save_path=None, mode='matrix'):
    """
    可视化一系列帧中密度矩阵的演化
    
    参数:
        frames: Frame对象列表或单个Frame对象
        system_name: 要可视化的系统名称（如果是None，则尝试找到第一个可用的系统）
        title: 可视化标题
        subsystems: 子系统维度，例如(15, 15)表示两个15维系统
        figsize: 图像大小
        save_path: 保存路径，如果为None则不保存
        mode: 可视化模式，'matrix'表示矩阵热图，'animation'表示动画
        
    返回:
        matplotlib.figure.Figure或matplotlib.animation.Animation: 根据mode返回相应对象
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import qutip as qt
    
    try:
        # 确保frames是列表
        if not isinstance(frames, list):
            frames = [frames]
        
        if not frames:
            print("没有提供帧数据")
            return None
        
        # 提取时间和密度矩阵
        times = []
        rho_evolution = []
        
        for frame in frames:
            times.append(frame.time)
            
            # 获取系统
            if system_name is None:
                # 尝试获取第一个可用的系统
                system_names = frame.get_all_system_names()
                if not system_names:
                    print(f"帧 {frame.frame_id} 中没有量子系统")
                    continue
                system_name = system_names[0]
            
            quantum_system = frame.get_quantum_system(system_name)
            if quantum_system is None:
                print(f"系统 {system_name} 在帧 {frame.frame_id} 中不存在")
                continue
            
            # 提取密度矩阵
            if isinstance(quantum_system, dict) and 'state' in quantum_system:
                state = quantum_system['state']
                rho_evolution.append(state)
            elif hasattr(quantum_system, 'state'):
                state = quantum_system.state
                rho_evolution.append(state)
            else:
                print(f"帧 {frame.frame_id} 中的系统 {system_name} 没有状态数据")
                print(f"量子系统类型: {type(quantum_system)}")
                if isinstance(quantum_system, dict):
                    print(f"键: {list(quantum_system.keys())}")
                continue
        
        if not rho_evolution:
            print("没有提取到密度矩阵数据")
            return None
        
        # 基于模式选择可视化方法
        if mode == 'matrix' and len(frames) == 1:
            # 单帧密度矩阵可视化 - 在独立窗口中显示
            plt.figure(figsize=figsize)
            fig = create_density_matrix_plot(
                rho_evolution[0], 
                title=f"{title} (t={times[0]:.3e})",
                subsystems=subsystems,
                figsize=figsize
            )
            
            if save_path:
                try:
                    fig.savefig(save_path, dpi=300, bbox_inches='tight')
                    print(f"图像已保存至 {save_path}")
                except Exception as e:
                    print(f"保存图像时出错: {e}")
            
            # 显示图形（无阻塞）
            plt.show(block=False)
            return fig
        
        elif mode == 'matrix' and len(frames) > 1:
            # 多帧比较可视化 (选择几个关键帧) - 在独立窗口中显示
            num_frames = min(4, len(frames))
            indices = np.linspace(0, len(frames)-1, num_frames, dtype=int)
            
            plt.figure(figsize=figsize)
            fig, axes = plt.subplots(2, num_frames, figsize=(figsize[0], figsize[1]*0.8))
            
            for i, idx in enumerate(indices):
                if hasattr(rho_evolution[idx], 'full'):
                    rho_data = rho_evolution[idx].full()
                else:
                    rho_data = np.array(rho_evolution[idx])
                    
                magnitude = np.abs(rho_data)
                phase = np.angle(rho_data)
                
                im_mag = axes[0, i].imshow(magnitude, cmap='viridis', interpolation='nearest')
                im_phase = axes[1, i].imshow(phase, cmap='hsv', interpolation='nearest', vmin=-np.pi, vmax=np.pi)
                
                axes[0, i].set_title(f't={times[idx]:.3e}')
                
                # 移除刻度标签，保持图像清晰
                for ax in [axes[0, i], axes[1, i]]:
                    ax.set_xticks([])
                    ax.set_yticks([])
            
            # 添加颜色条和标题
            plt.tight_layout()
            fig.colorbar(im_mag, ax=axes[0, :].ravel().tolist(), location='right', shrink=0.8)
            fig.colorbar(im_phase, ax=axes[1, :].ravel().tolist(), location='right', shrink=0.8)
            
            axes[0, 0].set_ylabel('幅值')
            axes[1, 0].set_ylabel('相位')
            
            fig.suptitle(title, fontsize=16)
            fig.subplots_adjust(top=0.9)
            
            if save_path:
                try:
                    fig.savefig(save_path, dpi=300, bbox_inches='tight')
                    print(f"图像已保存至 {save_path}")
                except Exception as e:
                    print(f"保存图像时出错: {e}")
            
            # 显示图形（无阻塞）
            plt.show(block=False)        
            return fig
        
        elif mode == 'animation':
            # 动画可视化 - 在独立窗口中显示
            ani = create_density_matrix_animation(
                rho_evolution,
                time_points=times,
                title=title,
                subsystems=subsystems,
                save_path=save_path
            )
            return ani
        
        else:
            print(f"不支持的可视化模式: {mode}")
            return None
    
    except Exception as e:
        import traceback
        print(f"可视化密度矩阵演化时出错: {e}")
        traceback.print_exc()
        return None

def verify_frame_files(directory_path, report_file=None, load_frames=False, plot_sample=False):
    """
    验证指定目录中的帧数据文件
    
    参数:
        directory_path (str): 包含帧数据文件的目录路径
        report_file (str, optional): 报告文件路径，若指定则将验证结果写入该文件
        load_frames (bool): 是否尝试加载帧数据
        plot_sample (bool): 是否绘制样本图像
        
    返回:
        dict: 包含验证结果和统计信息的字典
    """
    import os
    import json
    import numpy as np
    import glob
    import time
    from datetime import datetime
    
    # 初始化结果字典
    result = {
        "success": True,
        "total_files": 0,
        "valid_files": 0,
        "invalid_files": 0,
        "frames_count": 0,
        "files": [],
        "errors": [],
        "start_time": time.time()
    }
    
    # 检查目录是否存在
    if not os.path.exists(directory_path):
        result["success"] = False
        result["error"] = f"目录不存在: {directory_path}"
        return result
    
    # 查找所有帧数据文件
    frame_files = []
    
    # 检查frames目录
    frames_dir = os.path.join(directory_path, "frames")
    if os.path.exists(frames_dir):
        frame_files.extend(glob.glob(os.path.join(frames_dir, "frame_*.json")))
        frame_files.extend(glob.glob(os.path.join(frames_dir, "frame_*.npy")))
    
    # 检查根目录
    frame_files.extend(glob.glob(os.path.join(directory_path, "frame_*.json")))
    frame_files.extend(glob.glob(os.path.join(directory_path, "frame_*.npy")))
    
    # 统计文件数
    result["total_files"] = len(frame_files)
    
    if result["total_files"] == 0:
        result["success"] = False
        result["error"] = "未找到帧数据文件"
        return result
    
    # 记录开始验证
    print(f"找到 {result['total_files']} 个帧数据文件，开始验证...")
    
    # 分析每个文件
    for file_path in frame_files:
        file_info = {
            "path": file_path,
            "filename": os.path.basename(file_path),
            "size": os.path.getsize(file_path),
            "valid": False,
            "frames": 0,
            "format": os.path.splitext(file_path)[1],
            "timestamp": os.path.getmtime(file_path)
        }
        
        try:
            # 根据文件类型进行验证
            if file_path.endswith('.json'):
                # JSON格式文件
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                # 基本结构验证
                if not isinstance(data, dict):
                    raise ValueError("JSON文件格式错误，应为字典类型")
                
                if "frames" not in data:
                    raise ValueError("JSON文件缺少'frames'字段")
                
                frames = data["frames"]
                if not isinstance(frames, list):
                    raise ValueError("'frames'字段应为列表类型")
                
                file_info["frames"] = len(frames)
                
                # 验证每个帧的结构
                if load_frames and frames:
                    # 抽样检查，最多检查10个帧
                    sample_indices = np.linspace(0, len(frames)-1, min(10, len(frames))).astype(int)
                    
                    for idx in sample_indices:
                        frame = frames[idx]
                        if not isinstance(frame, dict):
                            raise ValueError(f"帧 {idx} 不是字典类型")
                        
                        # 检查必要字段
                        required_fields = ["t", "rho"]
                        for field in required_fields:
                            if field not in frame:
                                raise ValueError(f"帧 {idx} 缺少必要字段: {field}")
                        
                        # 检查rho字段是否为数组或嵌套列表
                        if not isinstance(frame["rho"], (list, dict)):
                            raise ValueError(f"帧 {idx} 的rho字段格式错误")
                        
                        # 如果是嵌套列表，检查其维度
                        if isinstance(frame["rho"], list) and frame["rho"]:
                            if isinstance(frame["rho"][0], list):
                                rows = len(frame["rho"])
                                if rows > 0:
                                    cols = len(frame["rho"][0])
                                    file_info["matrix_dim"] = f"{rows}x{cols}"
                
                file_info["valid"] = True
                result["valid_files"] += 1
                result["frames_count"] += file_info["frames"]
                
            elif file_path.endswith('.npy'):
                # NumPy格式文件
                data = np.load(file_path, allow_pickle=True)
                
                # 基本类型验证
                if not isinstance(data, (np.ndarray, dict)):
                    raise ValueError("NPY文件格式错误，应为数组或字典类型")
                
                # 如果是字典，检查frames字段
                if isinstance(data, dict):
                    if "frames" not in data:
                        raise ValueError("NPY文件缺少'frames'字段")
                    
                    frames = data["frames"]
                    if not isinstance(frames, (list, np.ndarray)):
                        raise ValueError("'frames'字段应为列表或数组类型")
                    
                    file_info["frames"] = len(frames)
                    
                    # 验证帧的结构
                    if load_frames and len(frames) > 0:
                        # 抽样检查
                        sample_indices = np.linspace(0, len(frames)-1, min(10, len(frames))).astype(int)
                        
                        for idx in sample_indices:
                            frame = frames[idx]
                            if not isinstance(frame, (dict, np.ndarray)):
                                raise ValueError(f"帧 {idx} 不是字典或数组类型")
                            
                            # 如果是字典，检查必要字段
                            if isinstance(frame, dict):
                                required_fields = ["t", "rho"]
                                for field in required_fields:
                                    if field not in frame:
                                        raise ValueError(f"帧 {idx} 缺少必要字段: {field}")
                                
                                # 检查rho字段
                                if "rho" in frame:
                                    rho = frame["rho"]
                                    if isinstance(rho, np.ndarray):
                                        file_info["matrix_dim"] = f"{rho.shape[0]}x{rho.shape[1]}"
                                        
                else:
                    # 如果是数组，假设它是帧的数组
                    file_info["frames"] = len(data)
                    
                    if load_frames and len(data) > 0:
                        # 抽样检查
                        sample_indices = np.linspace(0, len(data)-1, min(10, len(data))).astype(int)
                        
                        for idx in sample_indices:
                            frame = data[idx]
                            if not isinstance(frame, (dict, np.ndarray)):
                                raise ValueError(f"帧 {idx} 不是字典或数组类型")
                
                file_info["valid"] = True
                result["valid_files"] += 1
                result["frames_count"] += file_info["frames"]
                
            else:
                # 不支持的文件类型
                raise ValueError(f"不支持的文件类型: {os.path.splitext(file_path)[1]}")
            
        except Exception as e:
            # 记录错误
            file_info["error"] = str(e)
            result["invalid_files"] += 1
            result["errors"].append({
                "file": file_path,
                "error": str(e)
            })
            
            print(f"验证失败: {file_path}")
            print(f"  错误: {e}")
        
        # 添加文件信息
        result["files"].append(file_info)
    
    # 计算验证用时
    result["elapsed_time"] = time.time() - result["start_time"]
    
    # 输出摘要
    print("\n验证摘要:")
    print(f"总文件数: {result['total_files']}")
    print(f"有效文件: {result['valid_files']}")
    print(f"无效文件: {result['invalid_files']}")
    print(f"总帧数: {result['frames_count']}")
    print(f"验证用时: {result['elapsed_time']:.2f} 秒")
    
    # 如果需要，生成报告
    if report_file:
        try:
            # 创建目录（如果不存在）
            os.makedirs(os.path.dirname(report_file), exist_ok=True)
            
            with open(report_file, 'w') as f:
                f.write("=" * 50 + "\n")
                f.write(f"帧数据验证报告\n")
                f.write("=" * 50 + "\n\n")
                
                f.write(f"验证时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"目录: {directory_path}\n\n")
                
                f.write("摘要:\n")
                f.write(f"- 总文件数: {result['total_files']}\n")
                f.write(f"- 有效文件: {result['valid_files']}\n")
                f.write(f"- 无效文件: {result['invalid_files']}\n")
                f.write(f"- 总帧数: {result['frames_count']}\n")
                f.write(f"- 验证用时: {result['elapsed_time']:.2f} 秒\n\n")
                
                # 按有效性排序
                valid_files = [f for f in result["files"] if f["valid"]]
                invalid_files = [f for f in result["files"] if not f["valid"]]
                
                # 有效文件
                if valid_files:
                    f.write("有效文件:\n")
                    for file_info in valid_files:
                        f.write(f"- {file_info['filename']}:\n")
                        f.write(f"  * 格式: {file_info['format']}\n")
                        f.write(f"  * 大小: {file_info['size']} 字节\n")
                        f.write(f"  * 帧数: {file_info['frames']}\n")
                        if "matrix_dim" in file_info:
                            f.write(f"  * 矩阵维度: {file_info['matrix_dim']}\n")
                        f.write(f"  * 修改时间: {datetime.fromtimestamp(file_info['timestamp']).strftime('%Y-%m-%d %H:%M:%S')}\n")
                        f.write("\n")
                
                # 无效文件
                if invalid_files:
                    f.write("无效文件:\n")
                    for file_info in invalid_files:
                        f.write(f"- {file_info['filename']}:\n")
                        f.write(f"  * 错误: {file_info.get('error', '未知错误')}\n")
                        f.write("\n")
                
                # 错误详情
                if result["errors"]:
                    f.write("错误详情:\n")
                    for error in result["errors"]:
                        f.write(f"- {os.path.basename(error['file'])}:\n")
                        f.write(f"  * {error['error']}\n")
                        f.write("\n")
            
            print(f"报告已生成: {report_file}")
            result["report_file"] = report_file
        
        except Exception as e:
            print(f"报告生成失败: {e}")
            result["report_error"] = str(e)
    
    # 如果需要，绘制样本图像
    if plot_sample and result["valid_files"] > 0:
        try:
            import matplotlib.pyplot as plt
            from matplotlib.colors import ListedColormap
            import numpy as np
            
            print("绘制样本图像...")
            
            # 找到第一个有效文件
            valid_file = next((f for f in result["files"] if f["valid"]), None)
            
            if valid_file:
                # 根据文件类型加载数据
                if valid_file["path"].endswith('.json'):
                    with open(valid_file["path"], 'r') as f:
                        data = json.load(f)
                    frames = data["frames"]
                elif valid_file["path"].endswith('.npy'):
                    data = np.load(valid_file["path"], allow_pickle=True)
                    if isinstance(data, dict):
                        frames = data["frames"]
                    else:
                        frames = data
                
                # 确保有帧数据
                if frames and len(frames) > 0:
                    # 选择第一个帧
                    frame = frames[0]
                    
                    # 确保帧有正确的结构
                    if isinstance(frame, dict) and "rho" in frame:
                        rho = frame["rho"]
                        
                        # 如果rho是列表或数组，将其转换为数组
                        if isinstance(rho, list):
                            rho = np.array(rho)
                        elif isinstance(rho, dict):
                            # 如果是字典格式，尝试提取数据
                            if "data" in rho:
                                rho = np.array(rho["data"])
                            else:
                                print("无法识别的rho格式")
                                return result
                        
                        # 绘制密度矩阵
                        plt.figure(figsize=(10, 8))
                        
                        # 实部
                        plt.subplot(1, 2, 1)
                        plt.imshow(np.real(rho), cmap='RdBu_r', 
                                  vmin=-np.max(np.abs(np.real(rho))), 
                                  vmax=np.max(np.abs(np.real(rho))))
                        plt.colorbar(label='Real')
                        plt.title(f'Real part of ρ (t={frame.get("t", 0)})')
                        
                        # 虚部
                        plt.subplot(1, 2, 2)
                        plt.imshow(np.imag(rho), cmap='RdBu_r',
                                  vmin=-np.max(np.abs(np.imag(rho))), 
                                  vmax=np.max(np.abs(np.imag(rho))))
                        plt.colorbar(label='Imaginary')
                        plt.title(f'Imaginary part of ρ (t={frame.get("t", 0)})')
                        
                        plt.suptitle(f'验证样本: {os.path.basename(valid_file["path"])}')
                        plt.tight_layout()
                        
                        # 创建图像目录
                        if report_file:
                            img_dir = os.path.dirname(report_file)
                            img_path = os.path.join(img_dir, "sample_frame.png")
                            plt.savefig(img_path)
                            result["sample_image"] = img_path
                            print(f"样本图像已保存: {img_path}")
                        
                        plt.close()
        
        except Exception as e:
            print(f"样本图像绘制失败: {e}")
            result["plot_error"] = str(e)
    
    return result

def load_and_visualize_frames(file_path, frame_indices=None, plot_type="matrix", 
                             save_path=None, show=False, figsize=(15, 12)):
    """
    加载并可视化帧数据
    
    参数:
        file_path (str): 帧数据文件路径
        frame_indices (list): 要可视化的帧索引列表，为None则使用所有帧
        plot_type (str): 绘图类型，可为"matrix"(密度矩阵)或"wigner"(Wigner函数)
        save_path (str): 保存图像的路径，为None则不保存
        show (bool): 是否显示图像
        figsize (tuple): 图像大小
        
    返回:
        list: 加载的帧列表
    """
    import os
    import numpy as np
    import json
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    
    # 加载帧数据
    frames = []
    
    try:
        if file_path.endswith('.json'):
            with open(file_path, 'r') as f:
                data = json.load(f)
            frames = data.get("frames", [])
        elif file_path.endswith('.npy'):
            data = np.load(file_path, allow_pickle=True)
            if isinstance(data, dict):
                frames = data.get("frames", [])
            else:
                frames = data
        else:
            print(f"不支持的文件类型: {os.path.splitext(file_path)[1]}")
            return []
        
        # 如果没有指定帧索引，则使用所有帧
        if frame_indices is None:
            frame_indices = list(range(len(frames)))
        else:
            # 确保索引在有效范围内
            frame_indices = [i for i in frame_indices if 0 <= i < len(frames)]
        
        # 确保有有效的帧索引
        if not frame_indices:
            print("没有有效的帧索引")
            return frames
        
        # 创建图像
        n_frames = len(frame_indices)
        n_cols = min(3, n_frames)  # 每行最多3个图
        n_rows = (n_frames + n_cols - 1) // n_cols  # 向上取整
        
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(n_rows, n_cols, figure=fig)
        
        # 为每个选定的帧创建子图
        for i, idx in enumerate(frame_indices):
            row = i // n_cols
            col = i % n_cols
            
            frame = frames[idx]
            
            # 提取时间和密度矩阵
            if isinstance(frame, dict):
                t = frame.get("t", idx)
                rho = frame.get("rho", None)
            else:
                # 如果帧不是字典，则可能是直接的密度矩阵
                t = idx
                rho = frame
            
            # 转换rho为numpy数组（如果不是）
            if rho is not None:
                if isinstance(rho, list):
                    rho = np.array(rho)
                elif isinstance(rho, dict) and "data" in rho:
                    rho = np.array(rho["data"])
            
            # 跳过无效的rho
            if rho is None or not isinstance(rho, np.ndarray):
                plt.close()
        
        return loaded_frames
        
    except Exception as e:
        import traceback
        print(f"加载或可视化帧数据出错: {e}")
        traceback.print_exc()
        return None 

def load_simulation_frames(directory=None, simulation_id=None, limit=None):
    """
    从指定目录加载仿真帧数据
    
    参数:
        directory (str, optional): 帧数据目录，默认为results/visualization
        simulation_id (str, optional): 仿真ID，如果为None则加载最新的仿真
        limit (int, optional): 帧数量限制，如果为None则加载所有帧
        
    返回:
        list: 帧对象列表
    """
    import os
    import glob
    import h5py
    import numpy as np
    from core.quantum_state import Frame
    import time
    
    print("开始加载仿真帧数据...")
    
    # 确定目录
    if directory is None:
        directory = os.path.join("results", "visualization")
    
    if not os.path.exists(directory):
        print(f"目录不存在: {directory}")
        return []
    
    # 查找帧文件
    if simulation_id is not None:
        # 根据指定的仿真ID查找
        frame_files = glob.glob(os.path.join(directory, f"{simulation_id}_frames_*.h5"))
    else:
        # 查找所有帧文件
        frame_files = glob.glob(os.path.join(directory, "*_frames_*.h5"))
    
    if not frame_files:
        print(f"在目录 {directory} 中未找到帧文件")
        return []
    
    # 如果没有指定仿真ID，按修改时间排序并选择最新的仿真
    if simulation_id is None:
        frame_files = sorted(frame_files, key=os.path.getmtime, reverse=True)
        simulation_id = os.path.basename(frame_files[0]).split("_frames_")[0]
        print(f"选择最新的仿真: {simulation_id}")
    
    # 按时间顺序排序帧文件
    frame_parts = [(f, int(os.path.basename(f).split("_frames_")[1].split(".")[0])) 
                   for f in frame_files if simulation_id in f]
    frame_parts.sort(key=lambda x: x[1])
    
    # 加载帧
    frames = []
    start_time = time.time()
    
    try:
        for file_path, part_number in frame_parts:
            with h5py.File(file_path, 'r') as f:
                # 获取帧组
                if 'frames' not in f:
                    print(f"文件 {file_path} 中没有帧数据")
                    continue
                
                frames_group = f['frames']
                
                # 获取帧ID列表
                frame_ids = list(frames_group.keys())
                
                # 对帧ID按时间排序
                frame_times = []
                for frame_id in frame_ids:
                    frame_group = frames_group[frame_id]
                    if 'time' in frame_group.attrs:
                        frame_times.append((frame_id, frame_group.attrs['time']))
                    else:
                        # 从帧ID中提取时间戳（如果可能）
                        try:
                            # 假设帧ID格式为 frame_{timestamp}_{uuid}
                            timestamp = float(frame_id.split('_')[1])
                            frame_times.append((frame_id, timestamp))
                        except (IndexError, ValueError):
                            # 如果无法提取时间戳，则使用索引作为时间
                            frame_times.append((frame_id, float(len(frame_times))))
                
                # 按时间排序
                frame_times.sort(key=lambda x: x[1])
                
                # 限制帧数量
                if limit is not None and len(frames) + len(frame_times) > limit:
                    frame_times = frame_times[:limit - len(frames)]
                
                # 加载帧
                for frame_id, frame_time in frame_times:
                    frame_group = frames_group[frame_id]
                    
                    # 创建新帧
                    frame = Frame(frame_time)
                    
                    # 添加系统
                    if 'systems' in frame_group:
                        systems_group = frame_group['systems']
                        
                        for system_name in systems_group:
                            system_group = systems_group[system_name]
                            
                            # 提取状态数据
                            if 'state' in system_group:
                                state_dataset = system_group['state']
                                
                                # 提取实部和虚部
                                if 'real' in state_dataset and 'imag' in state_dataset:
                                    real_part = state_dataset['real'][()]
                                    imag_part = state_dataset['imag'][()]
                                    
                                    # 重构复数矩阵
                                    state = real_part + 1j * imag_part
                                else:
                                    # 如果没有分开存储，则直接读取
                                    state = state_dataset[()]
                                    
                                # 创建系统字典
                                system_data = {
                                    'system_name': system_name,
                                    'state': state
                                }
                                
                                # 添加额外属性
                                for attr_name, attr_value in system_group.attrs.items():
                                    system_data[attr_name] = attr_value
                                
                                # 添加到帧
                                frame.add_quantum_system(system_data)
                    
                    # 添加元数据
                    for attr_name, attr_value in frame_group.attrs.items():
                        frame.metadata[attr_name] = attr_value
                    
                    # 添加到帧列表
                    frames.append(frame)
                    
                    # 如果达到限制，停止加载
                    if limit is not None and len(frames) >= limit:
                        break
                
                # 如果达到限制，停止加载
                if limit is not None and len(frames) >= limit:
                    break
        
        print(f"成功加载 {len(frames)} 帧数据，用时 {time.time() - start_time:.2f} 秒")
        
        # 按时间排序
        frames.sort(key=lambda f: f.time)
        
        return frames
    
    except Exception as e:
        print(f"加载帧数据时出错: {e}")
        import traceback
        traceback.print_exc()
        return frames