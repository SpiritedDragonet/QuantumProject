import tkinter as tk
from tkinter import ttk, messagebox
import os
import sys
import threading
import time
import numpy as np
import queue
import platform

# 在所有matplotlib导入之前设置后端，避免创建多余窗口
import matplotlib
matplotlib.use('TkAgg')  # 强制使用TkAgg后端，防止创建独立窗口
# 关闭matplotlib的交互模式，防止创建多余窗口
matplotlib.interactive(False)
# 确保plt.show()不会启动事件循环
matplotlib.rcParams['figure.raise_window'] = False

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# 添加项目根目录到系统路径
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from config import create_default_config
from utils.visualization import ClassicalVisualization

# 直接设置SimHei作为默认中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题

# 通用Matplotlib设置
plt.rcParams.update({
    'font.size': 12,
    'figure.figsize': [10, 8],
    'figure.dpi': 100,
    'figure.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.format': 'png',
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'lines.linewidth': 1.5,
    'lines.markersize': 6
})

class SimpleSimulationData:
    """
    简单的仿真数据类，代替 SimulationData
    """
    def __init__(self):
        """初始化数据存储"""
        self.data = {}
        self.current_state = {
            'data_ready': False,
            'energy': 0.0
        }
    
    def update_data(self, key, value):
        """更新数据"""
        self.data[key] = value
    
    def get_data(self, key, default=None):
        """获取数据"""
        return self.data.get(key, default)
    
    def get_current_state(self):
        """获取当前状态"""
        return self.current_state
    
    def set_current_state(self, state):
        """设置当前状态"""
        self.current_state = state

class GeneralUI:
    """
    主界面UI类，作为整个系统的入口点
    """
    def __init__(self, root=None):
        """初始化主界面"""
        print("GeneralUI: Initializing...")
        self.simulation_thread = None
        self.stop_simulation = False
        self.sim_window = None  # 跟踪当前创建的仿真窗口
        
        # 如果没有提供根窗口，创建一个
        if root is None:
            print("GeneralUI: Creating Tkinter root window")
            self.root = tk.Tk()
            self.root.title("量子模拟系统")
            self.root.geometry("800x600")
            self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        else:
            self.root = root
        
        # 创建主框架
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)
        
        # 创建标题标签
        title_label = ttk.Label(self.main_frame, text="量子模拟系统", font=("SimHei", 24))
        title_label.pack(pady=20)
        
        # 创建描述标签
        description = """
        本系统用于模拟中性原子量子计算系统中的物理过程。
        
        点击下方按钮开始模拟。
        """
        desc_label = ttk.Label(self.main_frame, text=description, font=("SimHei", 12), wraplength=600, justify="center")
        desc_label.pack(pady=20)
        
        # 创建按钮框架
        button_frame = ttk.Frame(self.main_frame)
        button_frame.pack(pady=20)
        
        # 创建开始模拟按钮
        self.start_button = ttk.Button(button_frame, text="开始模拟", command=self.start_simulation)
        self.start_button.pack(side=tk.LEFT, padx=10)
        
        # 创建退出按钮
        exit_button = ttk.Button(button_frame, text="退出", command=self.on_close)
        exit_button.pack(side=tk.LEFT, padx=10)
        
        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # 模拟器实例
        self.simulator = None
        self.sim_ui = None
        self.sim_data = None
        
        # 帧数据结构初始化
        self.frame_data = {
            'frames': [],
            'current_index': 0,
            'total_frames': 0,
            'loaded': False,
            'source_file': None
        }
        
        print("GeneralUI: Initialization complete")
    
    def start_simulation(self):
        """启动量子仿真过程"""
        try:
            if self._simulate_boot_checks():
                self.status_var.set("仿真启动成功")
            else:
                self.status_var.set("系统自检失败，无法启动仿真")
        except Exception as e:
            print(f"启动仿真时出错: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("错误", f"启动仿真时出错: {e}")

    def _simulate_boot_checks(self):
        """
        检查系统环境和依赖项
        
        此方法在启动UI之前执行，检查Python版本、所需的库版本、
        目录结构以及关键文件是否存在。同时执行一次调试帧生成和验证，
        以确保系统正常运行。
        
        返回:
            bool: 检查是否通过
        """
        try:
            print("正在执行系统环境检查...")
            
            # 检查Python版本
            import sys
            python_version = sys.version.split()[0]
            print(f"Python版本: {python_version}")
            if not python_version.startswith('3'):
                print(f"错误: Python版本 {python_version} 不兼容，需要Python 3")
                return False
            
            # 检查必要的库
            try:
                import numpy as np
                numpy_version = np.__version__
                print(f"NumPy版本: {numpy_version}")
                
                import scipy
                scipy_version = scipy.__version__
                print(f"SciPy版本: {scipy_version}")
                
        import matplotlib
                matplotlib_version = matplotlib.__version__
                print(f"Matplotlib版本: {matplotlib_version}")
            except ImportError as e:
                print(f"错误: 缺少必要的库: {e}")
                return False
            
            # 检查目录结构
            import os
            # 检查NumPy版本
            numpy_version = np.__version__
            if not numpy_version:
                print("错误: 无法识别NumPy版本")
                return False
                
            # 检查SciPy版本
            scipy_version = scipy.__version__
            if not scipy_version:
                print("错误: 无法识别SciPy版本")
                return False
                
            # 检查Matplotlib版本
            matplotlib_version = matplotlib.__version__
            if not matplotlib_version:
                print("错误: 无法识别Matplotlib版本")
                return False
            
            # 检查目录结构
            import os
            project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            required_dirs = ['core', 'physics', 'protocols', 'utils', 'results']
            
            for dir_name in required_dirs:
                dir_path = os.path.join(project_root, dir_name)
                if not os.path.isdir(dir_path):
                    print(f"错误: 找不到必要目录 {dir_name}")
                    return False
            
            # 检查关键文件
            key_files = [
                os.path.join('utils', 'run_simulation.py'),
                os.path.join('core', 'quantum_state.py'),
                os.path.join('core', 'lindblad_solver.py')
            ]
            
            for file_path in key_files:
                abs_path = os.path.join(project_root, file_path)
                if not os.path.isfile(abs_path):
                    print(f"错误: 找不到关键文件 {file_path}")
            return False
    
            # 自检部分：生成调试帧和测试帧功能
            print("\n" + "="*60)
            print("系统自检：生成和验证调试帧")
            print("="*60)
            
            # 导入必要模块
            import time
            from datetime import datetime
            from config import create_default_config
            from utils.run_simulation import run_debug_simulation_with_frames
            from utils.visualization import verify_frame_files
            
            # 获取配置
            config = create_default_config()
            
            # 生成少量帧以进行快速测试（与原来的2000帧相比减少数量）
            test_frames_count = 100  # 自检时使用较少的帧数以提高速度
            test_duration = 100.0    # 自检时使用较短的时间以提高速度
            
            # 调整配置
            config['simulation_settings'] = config.get('simulation_settings', {})
            config['simulation_settings']['frames_count'] = test_frames_count
            config['simulation_settings']['total_time_ns'] = test_duration
            
            print(f"步骤1：生成 {test_frames_count} 帧仿真数据 (时长: {test_duration} ns)")
            
            start_time = time.time()
            
            # 运行调试模式仿真
            print(f"开始运行仿真自检...")
            simulation_dir = run_debug_simulation_with_frames(
                config,
                frames_count=test_frames_count,
                duration_ns=test_duration,
                is_self_check=True
            )
            
            elapsed = time.time() - start_time
            print(f"自检仿真完成！耗时 {elapsed:.2f} 秒")
            print(f"帧数据保存在: {simulation_dir}")
            
            # 步骤2：验证帧数据
            print("\n步骤2：验证帧数据文件")
            
            # 创建报告目录
            report_dir = os.path.join(simulation_dir, "reports")
            os.makedirs(report_dir, exist_ok=True)
            
            # 创建报告文件
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            report_file = os.path.join(report_dir, f"frame_verification_{timestamp}.txt")
            
            # 验证帧数据文件，但不绘图以提高速度
            verification_results = verify_frame_files(
                directory_path=simulation_dir,
                report_file=report_file,
                load_frames=True,
                plot_sample=False
            )
            
            # 检查验证结果
            if not verification_results["success"]:
                print(f"自检验证失败: {verification_results.get('error', '未知错误')}")
                return False
            else:
                print(f"自检验证成功! 找到 {verification_results.get('valid_files', 0)} 个有效帧文件，包含 {verification_results.get('frames_count', 0)} 帧数据")
            
            print("\n" + "="*60)
            print("系统自检完成！")
            print("="*60)
            
            return True
            
        except Exception as e:
            print(f"系统自检过程中出错: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _display_frame(self, frame_index):
        """
        显示指定索引的帧数据
        
        参数:
            frame_index: 帧索引
        """
        if frame_index < 0 or frame_index >= len(self.frame_data['frames']):
            print(f"无效的帧索引: {frame_index}")
            return
            
        frame = self.frame_data['frames'][frame_index]
        print(f"显示帧 {frame_index+1}/{len(self.frame_data['frames'])}")
        # 这里添加帧显示逻辑

    def load_frames(self, file_path):
        """
        从文件加载帧数据
        
        参数:
            file_path: 帧数据文件路径
            
        返回:
            bool: 是否成功加载
        """
        try:
            # 加载文件逻辑
            print(f"加载帧数据: {file_path}")
            return True
        except Exception as e:
            print(f"加载帧数据出错: {e}")
            return False
    
    def run(self):
        """运行主循环"""
        self.root.mainloop()
        
    def on_close(self):
        """
        关闭窗口时的处理函数
        """
        # 检查是否有仿真线程在运行
        if hasattr(self, 'simulation_thread') and self.simulation_thread and self.simulation_thread.is_alive():
            # 设置停止标志
            self.stop_simulation = True
            # 等待线程完成
            self.simulation_thread.join(timeout=0.5)  # 最多等待0.5秒
        
        # 销毁窗口
        self.root.destroy()
        
        print("GeneralUI: 应用已关闭")

def run_visualization(config=None, data_queue=None, load_path=None):
    """
    运行可视化界面
    
    参数:
        config: 配置对象
        data_queue: 数据队列
        load_path: 加载结果的路径
    """
    # matplotlib后端已在文件顶部设置，无需重复设置
    
    if config is None:
        config = create_default_config()
    
    print("运行可视化界面...")
    
    # 创建根窗口
    root = tk.Tk()
    
    # 创建UI实例
    ui = GeneralUI(root=root)  # 确保使用同一个root窗口
    
    # 如果提供了数据队列，设置监听
    if data_queue is not None:
        print("设置数据队列监听...")
        # 这里可以添加队列处理逻辑
    
    # 运行主循环
    ui.run()
    
    # 确保所有matplotlib窗口都被关闭
    plt.close('all')
