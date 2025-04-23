"""
系统状态记录与管理模块

本模块负责管理量子系统的演化过程中各个时间点的状态（帧），并提供存储和加载这些帧的功能。
主要功能包括：
- 记录系统演化的时间帧
- 将记录的帧保存到文件
- 从文件加载演化帧

函数调用关系：
```mermaid
graph TD
    A[SystemState.__init__] --> B[初始化系统状态管理器]
    C[record_frame] --> D[记录状态帧]
    E[dump_to_file] --> F[保存帧到文件]
    G[load_from_file] --> H[从文件加载帧]
```

数学背景：
- 帧集合 F = {(ρ_i, t_i, params_i) | i = 0,1,...,N-1}
- 每帧包含：密度矩阵 ρ_i、时间点 t_i 和附加参数 params_i

实现说明：
- 支持多个量子系统的帧记录
- 使用HDF5格式高效存储大型密度矩阵
"""

import numpy as np
import h5py
import json
import os
import pickle
from typing import List, Dict, Tuple, Optional, Union, Callable, Any
import datetime
import warnings
from pathlib import Path
import uuid
import qutip as qt

# 尝试导入core.quantum_state，如果导入失败则使用相对导入
try:
    from core.quantum_state import get_quantum_state, list_quantum_states, get_state_info, create_frame_data
except ImportError:
    from quantum_state import get_quantum_state, list_quantum_states, get_state_info, create_frame_data

# 帧处理相关函数
def _load_quantum_system_from_data(system_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    从帧数据中加载量子系统，将NumPy数组的密度矩阵转换为QuTiP的Qobj对象
    
    参数:
        system_data (Dict[str, Any]): 包含系统信息的字典
        
    返回:
        Dict[str, Any]: 更新后的系统信息字典
    """
    system_info = system_data.copy()
    
    # 获取维度和类型信息
    dims = system_info.get('dims', [])
    
    # 如果存在QuTiP的dims信息，使用它
    qutip_dims = system_info.get('qutip_dims', None)
    if qutip_dims is None:
        # 如果没有QuTiP的dims信息，构造标准的dims
        qutip_dims = [dims, dims]
    
    # 获取密度矩阵数据
    rho_data = system_info.get('rho', None)
    if rho_data is not None and not isinstance(rho_data, qt.Qobj):
        # 转换为QuTiP的Qobj对象
        system_info['rho'] = qt.Qobj(rho_data, dims=qutip_dims, type='oper')
    
    return system_info

class SystemState:
    """
    系统状态管理类
    
    该类负责记录量子系统在演化过程中的状态帧序列，并提供存储和加载这些帧的功能。
    
    主要属性：
        frames (List[Dict]): 存储系统演化的帧列表，每帧包含时间点、密度矩阵和附加参数
        metadata (Dict): 关于该演化过程的元数据
        simulation_id (str): 仿真过程的唯一标识符
    """
    
    def __init__(self, description: str = "", quantum_systems: Optional[List[str]] = None):
        """
        初始化SystemState对象
        
        参数：
            description (str): 本次仿真的描述信息
            quantum_systems (List[str], optional): 要关联的量子系统名称列表，默认为None表示所有系统
        """
        self.frames = []  # 帧列表，每个帧是包含rho、time和params的字典
        self.metadata = {
            "creation_time": datetime.datetime.now().isoformat(),
            "description": description,
            "frame_count": 0,
            "simulation_id": str(uuid.uuid4()),
            "version": "1.0"
        }
        self.simulation_id = self.metadata["simulation_id"]
        self.quantum_systems = quantum_systems if quantum_systems else list_quantum_states()
        self.last_saved_path = None
    
    def record_frame(self, t: float, params: Optional[Dict[str, Any]] = None, 
                    systems_to_record: Optional[List[str]] = None) -> None:
        """
        记录系统在时间t的状态帧
        
        该函数将系统在时间点t的状态和附加参数记录为一个帧，并添加到帧列表中。
        
        参数：
            t (float): 记录帧的时间点
            params (Dict[str, Any], optional): 附加参数信息，如期望值、纠缠度等
            systems_to_record (List[str], optional): 要记录的量子系统名称列表，默认为None表示记录所有系统
        
        数学表示：
            将帧 (ρ_systems, t, params) 添加到帧集合 F 中
            
        示例：
            # 记录系统在时间t=1.5的状态，并记录纠缠度等附加参数
            record_frame(1.5, {
                "concurrence": 0.85,
                "fidelity": 0.95,
                "photon_detection_probability": 0.1
            })
        """
        # 确定要记录的量子系统
        if systems_to_record is None:
            systems_to_record = self.quantum_systems
        
        # 构建帧数据
        frame = {
            "time": t,
            "timestamp": datetime.datetime.now().isoformat(),
            "quantum_systems": []
        }
        
        # 添加每个量子系统的状态信息
        for system_name in systems_to_record:
            try:
                state_info = get_state_info(system_name)
                frame["quantum_systems"].append(state_info)
            except Exception as e:
                warnings.warn(f"无法获取量子系统 '{system_name}' 的状态信息: {str(e)}")
        
        # 添加附加参数
        if params is not None:
            frame["params"] = params.copy()
        
        # 添加到帧列表
        self.frames.append(frame)
        self.metadata["frame_count"] = len(self.frames)
    
    def dump_to_file(self, filename: str, format: str = "hdf5", overwrite: bool = False) -> str:
        """
        将系统状态帧序列保存到文件
        
        将记录的帧序列序列化并保存到指定文件。支持HDF5、JSON和pickle格式。
        HDF5格式最适合大型密度矩阵，pickle适合快速保存，JSON适合跨平台兼容性。
        
        参数：
            filename (str): 保存文件的路径，如果没有扩展名，将根据format添加
            format (str): 文件格式，可选值为"hdf5"、"json"或"pickle"，默认为"hdf5"
            overwrite (bool): 是否覆盖已存在的文件，默认为False
            
        返回：
            str: 保存文件的完整路径
        """
        # 确保目录存在
        directory = os.path.dirname(filename)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)
        
        # 根据格式添加扩展名
        if not os.path.splitext(filename)[1]:
            if format == "hdf5":
                filename += ".h5"
            elif format == "json":
                filename += ".json"
            elif format == "pickle":
                filename += ".pkl"
            else:
                raise ValueError(f"不支持的文件格式: {format}")
        
        # 检查文件是否已存在
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(f"文件 {filename} 已存在。请使用overwrite=True参数覆盖现有文件。")
        
        # 更新元数据
        self.metadata["save_time"] = datetime.datetime.now().isoformat()
        self.metadata["frame_count"] = len(self.frames)
        self.metadata["file_format"] = format
        self.metadata["quantum_systems"] = self.quantum_systems
        
        # 根据不同格式保存
        if format == "hdf5":
            self._save_to_hdf5(filename)
        elif format == "json":
            self._save_to_json(filename)
        elif format == "pickle":
            self._save_to_pickle(filename)
        else:
            raise ValueError(f"不支持的文件格式: {format}")
        
        self.last_saved_path = filename
        return filename
    
    def _save_to_hdf5(self, filename: str) -> None:
        """
        将帧数据保存为HDF5格式
        
        HDF5格式适合保存大型密度矩阵和复杂的数据结构。
        
        参数：
            filename (str): HDF5文件路径
        """
        with h5py.File(filename, 'w') as f:
            # 保存元数据
            meta_group = f.create_group('metadata')
            for key, value in self.metadata.items():
                if isinstance(value, (str, int, float, bool)) or value is None:
                    meta_group.attrs[key] = value
                else:
                    # 对于复杂类型，转换为JSON字符串
                    meta_group.attrs[key] = json.dumps(value)
            
            # 创建帧组
            frames_group = f.create_group('frames')
            
            # 保存每一帧
            for i, frame in enumerate(self.frames):
                frame_group = frames_group.create_group(f'frame_{i:06d}')
                
                # 保存时间
                frame_group.attrs['time'] = frame['time']
                
                # 保存时间戳
                if 'timestamp' in frame:
                    frame_group.attrs['timestamp'] = frame['timestamp']
                
                # 保存量子系统信息
                if 'quantum_systems' in frame and frame['quantum_systems']:
                    systems_group = frame_group.create_group('quantum_systems')
                    for j, system_info in enumerate(frame['quantum_systems']):
                        system_group = systems_group.create_group(f'system_{j:03d}')
                        
                        # 保存系统名称
                        system_group.attrs['name'] = system_info['name']
                        
                        # 保存维度信息
                        system_group.attrs['total_dim'] = system_info['total_dim']
                        system_group.create_dataset('dims', data=system_info['dims'])
                        system_group.create_dataset('types', data=system_info['types'])
                        
                        # 保存QuTiP的dims信息
                        if 'qutip_dims' in system_info:
                            system_group.attrs['qutip_dims'] = json.dumps(system_info['qutip_dims'])
                        
                        # 保存密度矩阵
                        # 如果是Qobj对象，转换为NumPy数组
                        if isinstance(system_info['rho'], qt.Qobj):
                            density_matrix = system_info['rho'].full()
                        else:
                            density_matrix = system_info['rho']
                        system_group.create_dataset('density_matrix', data=density_matrix)
                        
                        # 保存子系统结构描述
                        if 'subsystem_structure' in system_info:
                            system_group.attrs['subsystem_structure'] = json.dumps(system_info['subsystem_structure'])
                
                # 保存附加参数
                if 'params' in frame and frame['params']:
                    params_group = frame_group.create_group('params')
                    self._save_dict_to_hdf5_group(params_group, frame['params'])
    
    def _save_dict_to_hdf5_group(self, group, data_dict: Dict) -> None:
        """
        将字典保存到HDF5组
        
        递归地将字典中的数据保存到HDF5组。
        
        参数：
            group: HDF5组对象
            data_dict (Dict): 要保存的字典
        """
        for key, value in data_dict.items():
            if value is None:
                group.attrs[key] = 'None'
            elif isinstance(value, (str, int, float, bool)):
                group.attrs[key] = value
            elif isinstance(value, np.ndarray):
                group.create_dataset(key, data=value)
            elif isinstance(value, qt.Qobj):
                # 保存QuTiP的Qobj对象为NumPy数组
                group.create_dataset(key, data=value.full())
                group.attrs[f"{key}_qutip_dims"] = json.dumps(value.dims)
                group.attrs[f"{key}_qutip_type"] = value.type
            elif isinstance(value, dict):
                subgroup = group.create_group(key)
                self._save_dict_to_hdf5_group(subgroup, value)
            elif isinstance(value, list):
                if value and all(isinstance(x, (int, float, str, bool)) for x in value):
                    # 简单类型的列表
                    try:
                        # 尝试转换为同质数组
                        array_value = np.array(value)
                        group.create_dataset(key, data=array_value)
                    except:
                        # 异质列表转为JSON
                        group.attrs[key] = json.dumps(value)
                else:
                    # 复杂类型的列表转为JSON
                    group.attrs[key] = json.dumps(value)
            else:
                # 其他类型转为字符串
                group.attrs[key] = str(value)
    
    def _save_to_json(self, filename: str) -> None:
        """
        将帧数据保存为JSON格式
        
        JSON格式适合保存小规模系统的数据，具有良好的可读性和跨平台兼容性。
        
        参数：
            filename (str): JSON文件路径
        """
        # 创建可序列化的数据结构
        data = {
            "metadata": self.metadata.copy(),
            "frames": []
        }
        
        # 处理每一帧
        for frame in self.frames:
            frame_data = {
                "time": frame["time"]
            }
            
            # 处理时间戳
            if "timestamp" in frame:
                frame_data["timestamp"] = frame["timestamp"]
            
            # 处理量子系统信息
            if "quantum_systems" in frame and frame["quantum_systems"]:
                frame_data["quantum_systems"] = []
                for system_info in frame["quantum_systems"]:
                    system_data = {
                        "name": system_info["name"],
                        "total_dim": system_info["total_dim"],
                        "dims": system_info["dims"],
                        "types": system_info["types"]
                    }
                    
                    # 处理密度矩阵
                    if isinstance(system_info["rho"], qt.Qobj):
                        density_matrix = system_info["rho"].full()
                        system_data["qutip_dims"] = system_info["rho"].dims
                    else:
                        density_matrix = system_info["rho"]
                        if "qutip_dims" in system_info:
                            system_data["qutip_dims"] = system_info["qutip_dims"]
                    
                    system_data["density_matrix"] = {
                        "real": density_matrix.real.tolist(),
                        "imag": density_matrix.imag.tolist(),
                        "shape": density_matrix.shape
                    }
                    
                    # 处理子系统结构描述
                    if "subsystem_structure" in system_info:
                        system_data["subsystem_structure"] = system_info["subsystem_structure"]
                    
                    frame_data["quantum_systems"].append(system_data)
            
            # 处理附加参数
            if "params" in frame and frame["params"]:
                frame_data["params"] = self._convert_to_json_serializable(frame["params"])
            
            data["frames"].append(frame_data)
        
        # 写入文件
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
    
    def _convert_to_json_serializable(self, obj):
        """
        将对象转换为JSON可序列化的格式
        
        递归地将包含numpy数组、复数等不可JSON序列化的对象转换为可序列化的格式。
        
        参数：
            obj: 要转换的对象
            
        返回：
            转换后的可JSON序列化对象
        """
        if obj is None:
            return None
        elif isinstance(obj, (int, float, str, bool)):
            return obj
        elif isinstance(obj, complex):
            return {"real": obj.real, "imag": obj.imag, "_type": "complex"}
        elif isinstance(obj, np.ndarray):
            if np.issubdtype(obj.dtype, np.complexfloating):
                return {
                    "real": obj.real.tolist(),
                    "imag": obj.imag.tolist(),
                    "shape": obj.shape,
                    "_type": "complex_array"
                }
            else:
                return {
                    "data": obj.tolist(),
                    "shape": obj.shape,
                    "_type": "array"
                }
        elif isinstance(obj, qt.Qobj):
            return {
                "real": obj.full().real.tolist(),
                "imag": obj.full().imag.tolist(),
                "shape": obj.shape,
                "dims": obj.dims,
                "type": obj.type,
                "_type": "qutip_qobj"
            }
        elif isinstance(obj, dict):
            return {k: self._convert_to_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_to_json_serializable(item) for item in obj]
        else:
            return str(obj)
    
    def _save_to_pickle(self, filename: str) -> None:
        """
        将帧数据保存为pickle格式
        
        pickle格式适合快速保存数据，但可能存在跨平台和版本兼容性问题。
        
        参数：
            filename (str): pickle文件路径
        """
        data = {
            "metadata": self.metadata,
            "frames": self.frames
        }
        with open(filename, 'wb') as f:
            pickle.dump(data, f)
    
    def load_from_file(self, filename: str) -> None:
        """
        从文件加载系统状态帧序列
        
        从先前保存的文件加载系统状态帧序列，支持HDF5、JSON和pickle格式。
        
        参数：
            filename (str): 要加载的文件路径
        """
        if not os.path.exists(filename):
            raise FileNotFoundError(f"文件 {filename} 不存在")
        
        # 根据文件扩展名确定格式
        ext = os.path.splitext(filename)[1].lower()
        if ext == '.h5' or ext == '.hdf5':
            self._load_from_hdf5(filename)
        elif ext == '.json':
            self._load_from_json(filename)
        elif ext == '.pkl' or ext == '.pickle':
            self._load_from_pickle(filename)
        else:
            raise ValueError(f"不支持的文件格式: {ext}")
        
        self.last_saved_path = filename
    
    def _load_from_hdf5(self, filename: str) -> None:
        """
        从HDF5文件加载帧数据
        
        参数：
            filename (str): HDF5文件路径
        """
        with h5py.File(filename, 'r') as f:
            # 加载元数据
            self.metadata = {}
            meta_group = f['metadata']
            for key in meta_group.attrs:
                value = meta_group.attrs[key]
                # 尝试解析JSON字符串
                if isinstance(value, str) and (value.startswith('{') or value.startswith('[')):
                    try:
                        self.metadata[key] = json.loads(value)
                    except:
                        self.metadata[key] = value
                else:
                    self.metadata[key] = value
            
            # 设置仿真ID和量子系统列表
            self.simulation_id = self.metadata.get("simulation_id", str(uuid.uuid4()))
            self.quantum_systems = self.metadata.get("quantum_systems", list_quantum_states())
            
            # 加载帧
            self.frames = []
            frames_group = f['frames']
            for frame_name in sorted(frames_group.keys()):
                frame_group = frames_group[frame_name]
                
                # 创建帧字典
                frame = {}
                
                # 加载时间
                frame['time'] = frame_group.attrs['time']
                
                # 加载时间戳
                if 'timestamp' in frame_group.attrs:
                    frame['timestamp'] = frame_group.attrs['timestamp']
                
                # 加载量子系统信息
                if 'quantum_systems' in frame_group:
                    systems_group = frame_group['quantum_systems']
                    frame['quantum_systems'] = []
                    
                    for system_name in sorted(systems_group.keys()):
                        system_group = systems_group[system_name]
                        
                        # 创建系统信息字典
                        system_info = {}
                        
                        # 加载系统名称
                        system_info['name'] = system_group.attrs['name']
                        
                        # 加载维度信息
                        system_info['total_dim'] = system_group.attrs['total_dim']
                        system_info['dims'] = system_group['dims'][()]
                        system_info['types'] = system_group['types'][()]
                        
                        # 加载QuTiP的dims信息
                        if 'qutip_dims' in system_group.attrs:
                            qutip_dims_str = system_group.attrs['qutip_dims']
                            try:
                                system_info['qutip_dims'] = json.loads(qutip_dims_str)
                            except:
                                pass
                        
                        # 加载密度矩阵
                        density_matrix = system_group['density_matrix'][()]
                        
                        # 转换为QuTiP的Qobj对象
                        if isinstance(system_info['rho'], qt.Qobj):
                            density_matrix = system_info['rho'].full()
                        else:
                            density_matrix = system_info['rho']
                        
                        system_info['rho'] = qt.Qobj(density_matrix, dims=system_info['qutip_dims'], type='oper')
                        
                        # 加载子系统结构描述
                        if 'subsystem_structure' in system_group.attrs:
                            subsystem_structure = system_group.attrs['subsystem_structure']
                            if isinstance(subsystem_structure, str):
                                try:
                                    system_info['subsystem_structure'] = json.loads(subsystem_structure)
                                except:
                                    system_info['subsystem_structure'] = subsystem_structure
                            else:
                                system_info['subsystem_structure'] = subsystem_structure
                        
                        frame['quantum_systems'].append(system_info)
                
                # 加载附加参数
                if 'params' in frame_group:
                    frame['params'] = self._load_dict_from_hdf5_group(frame_group['params'])
                
                self.frames.append(frame)
    
    def _load_dict_from_hdf5_group(self, group) -> Dict:
        """
        从HDF5组加载字典
        
        递归地从HDF5组加载数据到字典。
        
        参数：
            group: HDF5组对象
            
        返回：
            Dict: 加载的字典
        """
        result = {}
        
        # 加载属性
        for key in group.attrs:
            value = group.attrs[key]
            if value == 'None':
                result[key] = None
            elif isinstance(value, str) and (value.startswith('{') or value.startswith('[')):
                try:
                    result[key] = json.loads(value)
                except:
                    result[key] = value
            else:
                result[key] = value
        
        # 检查是否有QuTiP对象的元数据
        qutip_keys = set()
        for key in group.attrs:
            if key.endswith('_qutip_dims') or key.endswith('_qutip_type'):
                base_key = key.rsplit('_', 1)[0]
                qutip_keys.add(base_key)
        
        # 加载数据集和子组
        for key in group:
            item = group[key]
            if isinstance(item, h5py.Dataset):
                data = item[()]
                
                # 检查是否是QuTiP对象
                if key in qutip_keys:
                    dims_key = f"{key}_qutip_dims"
                    type_key = f"{key}_qutip_type"
                    if dims_key in group.attrs:
                        dims = json.loads(group.attrs[dims_key])
                        qobj_type = group.attrs.get(type_key, 'oper')
                        result[key] = qt.Qobj(data, dims=dims, type=qobj_type)
                    else:
                        result[key] = data
                else:
                    result[key] = data
            elif isinstance(item, h5py.Group):
                result[key] = self._load_dict_from_hdf5_group(item)
        
        return result
    
    def _load_from_json(self, filename: str) -> None:
        """
        从JSON文件加载帧数据
        
        参数：
            filename (str): JSON文件路径
        """
        with open(filename, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # 加载元数据
        self.metadata = data["metadata"]
        self.simulation_id = self.metadata.get("simulation_id", str(uuid.uuid4()))
        self.quantum_systems = self.metadata.get("quantum_systems", list_quantum_states())
        
        # 加载帧
        self.frames = []
        for frame_data in data["frames"]:
            frame = {
                "time": frame_data["time"]
            }
            
            # 加载时间戳
            if "timestamp" in frame_data:
                frame["timestamp"] = frame_data["timestamp"]
            
            # 加载量子系统信息
            if "quantum_systems" in frame_data and frame_data["quantum_systems"]:
                frame["quantum_systems"] = []
                
                for system_data in frame_data["quantum_systems"]:
                    system_info = {
                        "name": system_data["name"],
                        "total_dim": system_data["total_dim"],
                        "dims": system_data["dims"],
                        "types": system_data["types"]
                    }
                    
                    # 加载密度矩阵
                    dm_data = system_data["density_matrix"]
                    real_part = np.array(dm_data["real"])
                    imag_part = np.array(dm_data["imag"])
                    
                    # 加载QuTiP的dims信息
                    if "qutip_dims" in system_data:
                        system_info["qutip_dims"] = system_data["qutip_dims"]
                        qutip_dims = system_data["qutip_dims"]
                        # 创建QuTiP对象
                        system_info["rho"] = qt.Qobj(real_part + 1j * imag_part, 
                                                   dims=qutip_dims, 
                                                   type='oper')
                    else:
                        # 使用标准dims创建QuTiP对象
                        dims = system_info["dims"]
                        system_info["rho"] = qt.Qobj(real_part + 1j * imag_part, 
                                                   dims=[dims, dims], 
                                                   type='oper')
                    
                    # 加载子系统结构描述
                    if "subsystem_structure" in system_data:
                        system_info["subsystem_structure"] = system_data["subsystem_structure"]
                    
                    frame["quantum_systems"].append(system_info)
            
            # 加载附加参数
            if "params" in frame_data and frame_data["params"]:
                frame["params"] = self._convert_from_json_serializable(frame_data["params"])
            
            self.frames.append(frame)
    
    def _convert_from_json_serializable(self, obj):
        """
        将JSON反序列化的对象转换回原始类型
        
        递归地将反序列化的JSON对象转换回原始类型，如复数和numpy数组。
        
        参数：
            obj: 要转换的对象
            
        返回：
            转换后的对象
        """
        if obj is None:
            return None
        elif isinstance(obj, (int, float, str, bool)):
            return obj
        elif isinstance(obj, dict):
            if "_type" in obj:
                if obj["_type"] == "complex":
                    return complex(obj["real"], obj["imag"])
                elif obj["_type"] == "complex_array":
                    return np.array(obj["real"]) + 1j * np.array(obj["imag"])
                elif obj["_type"] == "array":
                    return np.array(obj["data"])
                elif obj["_type"] == "qutip_qobj":
                    data = np.array(obj["real"]) + 1j * np.array(obj["imag"])
                    return qt.Qobj(data, dims=obj.get("dims"), type=obj.get("type", "oper"))
            return {k: self._convert_from_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_from_json_serializable(item) for item in obj]
        else:
            return obj
    
    def _load_from_pickle(self, filename: str) -> None:
        """
        从pickle文件加载帧数据
        
        参数：
            filename (str): pickle文件路径
        """
        with open(filename, 'rb') as f:
            data = pickle.load(f)
        
        self.metadata = data["metadata"]
        self.simulation_id = self.metadata.get("simulation_id", str(uuid.uuid4()))
        self.quantum_systems = self.metadata.get("quantum_systems", list_quantum_states())
        
        # 加载帧并确保密度矩阵是QuTiP的Qobj对象
        self.frames = []
        for frame in data["frames"]:
            # 处理量子系统
            if "quantum_systems" in frame:
                for i, system_info in enumerate(frame["quantum_systems"]):
                    # 将NumPy数组转换为QuTiP的Qobj对象
                    frame["quantum_systems"][i] = _load_quantum_system_from_data(system_info)
            
            self.frames.append(frame)

# 全局单例系统状态对象
_system_state_instance = None

def get_system_state(description: str = "", quantum_systems: Optional[List[str]] = None) -> SystemState:
    """
    获取全局系统状态管理器实例
    
    参数：
        description (str): 系统状态描述
        quantum_systems (List[str], optional): 要关联的量子系统名称列表，默认为None表示所有系统
        
    返回：
        SystemState: 全局系统状态管理器实例
    """
    global _system_state_instance
    if _system_state_instance is None:
        _system_state_instance = SystemState(description, quantum_systems)
    return _system_state_instance

def record_frame(t: float, params: Optional[Dict[str, Any]] = None,
                systems_to_record: Optional[List[str]] = None) -> None:
    """
    记录系统在时间t的状态帧（全局函数）
    
    参数：
        t (float): 记录帧的时间点
        params (Dict[str, Any], optional): 附加参数信息
        systems_to_record (List[str], optional): 要记录的量子系统名称列表，默认为None表示记录所有系统
    """
    system_state = get_system_state()
    system_state.record_frame(t, params, systems_to_record)

def dump_to_file(filename: str, format: str = "hdf5", overwrite: bool = False) -> str:
    """
    将系统状态帧序列保存到文件（全局函数）
    
    参数：
        filename (str): 保存文件的路径
        format (str): 文件格式，可选值为"hdf5"、"json"或"pickle"，默认为"hdf5"
        overwrite (bool): 是否覆盖已存在的文件，默认为False
            
    返回：
        str: 保存文件的完整路径
    """
    system_state = get_system_state()
    return system_state.dump_to_file(filename, format, overwrite)

def load_from_file(filename: str) -> None:
    """
    从文件加载系统状态帧序列（全局函数）
    
    参数：
        filename (str): 要加载的文件路径
    """
    system_state = get_system_state()
    system_state.load_from_file(filename)

def create_frame_data(self, t: float, additional_data: Optional[Dict] = None) -> Dict:
    """
    创建包含当前仿真状态信息的帧数据
    
    构造一个包含当前时间点的量子系统状态和附加信息的帧数据字典。
    
    参数：
        t (float): 当前时间
        additional_data (Dict, optional): 任何附加的想要记录的数据
        
    返回：
        Dict: 包含完整帧数据的字典，包括：
            - 't': 时间点
            - 'quantum_states': 包含所有被记录的量子系统状态信息的列表
            - 其他任何通过additional_data提供的数据
    """
    frame = {
        't': t,
        'quantum_states': []
    }
    
    # 收集所有量子系统的信息
    for system_name in self.quantum_system_names:
        try:
            # 获取量子系统
            system = quantum_state.get_quantum_state(system_name)
            
            # 获取状态信息，包括密度矩阵、哈密顿量和Lindblad算符
            state_info = system.get_state_info()
            
            # 添加到帧中
            frame['quantum_states'].append(state_info)
        except Exception as e:
            warnings.warn(f"无法获取量子系统 '{system_name}' 的状态信息: {str(e)}")
    
    # 添加额外数据（如果有）
    if additional_data:
        for key, value in additional_data.items():
            frame[key] = value
    
    return frame

def _save_frame_to_hdf5(self, group, frame: Dict, frame_index: int) -> None:
    """
    将单个帧保存到HDF5文件的指定组中
    
    参数：
        group (h5py.Group): HDF5文件中的组对象
        frame (Dict): 帧数据字典
        frame_index (int): 帧索引
    """
    # 创建帧子组
    frame_group = group.create_group(f'frame_{frame_index}')
    
    # 保存时间
    frame_group.attrs['t'] = frame['t']
    
    # 保存量子态信息
    if 'quantum_states' in frame:
        states_group = frame_group.create_group('quantum_states')
        
        for i, state_info in enumerate(frame['quantum_states']):
            state_group = states_group.create_group(f'state_{i}')
            
            # 保存状态基本信息
            for key in ['name', 'total_dim']:
                if key in state_info:
                    state_group.attrs[key] = state_info[key]
            
            # 保存数组类型数据
            for key in ['dims', 'types']:
                if key in state_info:
                    state_group.create_dataset(key, data=state_info[key])
            
            # 保存密度矩阵
            if 'rho' in state_info:
                state_group.create_dataset('rho', data=state_info['rho'])
                # 保存QuTiP dims信息
                if 'qutip_dims' in state_info:
                    dims_group = state_group.create_group('qutip_dims')
                    for j, dim_list in enumerate(state_info['qutip_dims']):
                        dims_group.create_dataset(f'dim_{j}', data=dim_list)
            
            # 保存哈密顿量（如果存在）
            if 'hamiltonian' in state_info:
                state_group.create_dataset('hamiltonian', data=state_info['hamiltonian'])
                # 保存QuTiP dims信息
                if 'hamiltonian_qutip_dims' in state_info:
                    h_dims_group = state_group.create_group('hamiltonian_qutip_dims')
                    for j, dim_list in enumerate(state_info['hamiltonian_qutip_dims']):
                        h_dims_group.create_dataset(f'dim_{j}', data=dim_list)
            
            # 保存Lindblad算符（如果存在）
            if 'lindblad_ops' in state_info:
                lindblad_group = state_group.create_group('lindblad_ops')
                for j, op_data in enumerate(state_info['lindblad_ops']):
                    lindblad_group.create_dataset(f'op_{j}', data=op_data)
                
                # 保存QuTiP dims信息
                if 'lindblad_qutip_dims' in state_info:
                    lindblad_dims_group = state_group.create_group('lindblad_qutip_dims')
                    for j, op_dims in enumerate(state_info['lindblad_qutip_dims']):
                        op_dims_group = lindblad_dims_group.create_group(f'op_{j}_dims')
                        for k, dim_list in enumerate(op_dims):
                            op_dims_group.create_dataset(f'dim_{k}', data=dim_list)
            
            # 保存子系统结构描述
            if 'subsystem_structure' in state_info:
                state_group.attrs['subsystem_structure'] = json.dumps(state_info['subsystem_structure'])
    
    # 保存额外数据
    for key, value in frame.items():
        if key not in ['t', 'quantum_states']:
            # 如果是简单数据类型，保存为属性
            if isinstance(value, (int, float, str, bool)) or value is None:
                frame_group.attrs[key] = value
            # 如果是数组，保存为数据集
            elif isinstance(value, (list, np.ndarray)):
                # 对于不能直接保存的数据类型，尝试JSON序列化
                try:
                    if isinstance(value, list) and not isinstance(value[0], (int, float, bool, str)):
                        frame_group.attrs[key] = json.dumps(value)
                    else:
                        frame_group.create_dataset(key, data=np.array(value))
                except (TypeError, ValueError):
                    warnings.warn(f"无法保存帧额外数据项 '{key}'，数据类型不兼容")
            # 对于字典，转为JSON
            elif isinstance(value, dict):
                try:
                    frame_group.attrs[key] = json.dumps(value)
                except (TypeError, ValueError):
                    warnings.warn(f"无法保存帧额外数据项 '{key}'，无法JSON序列化")

def _load_frame_from_hdf5(self, group) -> Dict:
    """
    从HDF5文件的指定组中加载单个帧
    
    参数：
        group (h5py.Group): HDF5文件中的组对象
        
    返回：
        Dict: 加载的帧数据字典
    """
    frame = {}
    
    # 加载时间
    frame['t'] = group.attrs['t']
    
    # 加载量子态信息
    if 'quantum_states' in group:
        frame['quantum_states'] = []
        states_group = group['quantum_states']
        
        for state_key in sorted(states_group.keys()):
            state_group = states_group[state_key]
            state_info = {}
            
            # 加载状态基本信息
            for key in state_group.attrs:
                state_info[key] = state_group.attrs[key]
            
            # 加载数组类型数据
            for key in ['dims', 'types']:
                if key in state_group:
                    state_info[key] = state_group[key][()]
            
            # 加载密度矩阵
            if 'rho' in state_group:
                state_info['rho'] = state_group['rho'][()]
                
                # 加载QuTiP dims信息
                if 'qutip_dims' in state_group:
                    dims_group = state_group['qutip_dims']
                    qutip_dims = []
                    for i in range(len(dims_group)):
                        dim_key = f'dim_{i}'
                        if dim_key in dims_group:
                            qutip_dims.append(dims_group[dim_key][()])
                    state_info['qutip_dims'] = qutip_dims
            
            # 加载哈密顿量（如果存在）
            if 'hamiltonian' in state_group:
                state_info['hamiltonian'] = state_group['hamiltonian'][()]
                
                # 加载QuTiP dims信息
                if 'hamiltonian_qutip_dims' in state_group:
                    h_dims_group = state_group['hamiltonian_qutip_dims']
                    h_qutip_dims = []
                    for i in range(len(h_dims_group)):
                        dim_key = f'dim_{i}'
                        if dim_key in h_dims_group:
                            h_qutip_dims.append(h_dims_group[dim_key][()])
                    state_info['hamiltonian_qutip_dims'] = h_qutip_dims
            
            # 加载Lindblad算符（如果存在）
            if 'lindblad_ops' in state_group:
                lindblad_group = state_group['lindblad_ops']
                lindblad_ops = []
                
                # 确定算符数量
                op_count = sum(1 for key in lindblad_group.keys() if key.startswith('op_'))
                
                # 按顺序加载所有算符
                for i in range(op_count):
                    op_key = f'op_{i}'
                    if op_key in lindblad_group:
                        lindblad_ops.append(lindblad_group[op_key][()])
                
                state_info['lindblad_ops'] = lindblad_ops
                
                # 加载QuTiP dims信息
                if 'lindblad_qutip_dims' in state_group:
                    lindblad_dims_group = state_group['lindblad_qutip_dims']
                    lindblad_qutip_dims = []
                    
                    # 确定算符数量
                    op_dims_count = sum(1 for key in lindblad_dims_group.keys() if key.endswith('_dims'))
                    
                    # 按顺序加载所有算符的dims信息
                    for i in range(op_dims_count):
                        op_dims_key = f'op_{i}_dims'
                        if op_dims_key in lindblad_dims_group:
                            op_dims_group = lindblad_dims_group[op_dims_key]
                            op_dims = []
                            
                            # 加载该算符的dims列表
                            for j in range(len(op_dims_group)):
                                dim_key = f'dim_{j}'
                                if dim_key in op_dims_group:
                                    op_dims.append(op_dims_group[dim_key][()])
                            
                            lindblad_qutip_dims.append(op_dims)
                    
                    state_info['lindblad_qutip_dims'] = lindblad_qutip_dims
            
            # 加载子系统结构描述
            if 'subsystem_structure' in state_group.attrs:
                try:
                    state_info['subsystem_structure'] = json.loads(state_group.attrs['subsystem_structure'])
                except:
                    state_info['subsystem_structure'] = state_group.attrs['subsystem_structure']
            
            frame['quantum_states'].append(state_info)
    
    # 加载额外数据
    for key in group.attrs:
        if key != 't':
            value = group.attrs[key]
            # 尝试解析JSON字符串
            if isinstance(value, str) and (value.startswith('{') or value.startswith('[')):
                try:
                    value = json.loads(value)
                except json.JSONDecodeError:
                    pass
            frame[key] = value
    
    # 加载数据集形式的额外数据
    for key in group.keys():
        if key != 'quantum_states' and not key.startswith('frame_'):
            try:
                frame[key] = group[key][()]
            except:
                pass
    
    return frame

def restore_quantum_states_from_frame(self, frame_index: int) -> None:
    """
    从指定帧恢复量子系统状态
    
    从存储的帧数据中恢复所有量子系统的状态，包括密度矩阵、哈密顿量和Lindblad算符
    
    参数：
        frame_index (int): 要恢复的帧索引
    """
    if not self.frames or frame_index < 0 or frame_index >= len(self.frames):
        raise ValueError(f"无效的帧索引: {frame_index}，有效范围为0-{len(self.frames)-1}")
    
    frame = self.frames[frame_index]
    
    if 'quantum_states' not in frame:
        raise ValueError(f"帧 {frame_index} 不包含量子系统状态信息")
    
    # 恢复每个量子系统的状态
    for state_info in frame['quantum_states']:
        if 'name' not in state_info:
            warnings.warn("帧中的量子系统状态缺少名称，无法恢复")
            continue
        
        system_name = state_info['name']
        
        # 检查该量子系统是否已存在，如果不存在则创建
        try:
            system = quantum_state.get_quantum_state(system_name)
        except ValueError:
            system = quantum_state.register_quantum_state(quantum_state.QuantumState(system_name))
        
        # 初始化希尔伯特空间
        if 'dims' in state_info and 'types' in state_info:
            # 转换类型数值为枚举
            types = [quantum_state.SubsystemType(t) for t in state_info['types']]
            system.init_hilbert_space(state_info['dims'], types)
        
        # 恢复密度矩阵
        if 'rho' in state_info:
            rho_data = state_info['rho']
            
            # 如果有QuTiP dims信息，使用它创建Qobj
            if 'qutip_dims' in state_info:
                rho = qt.Qobj(rho_data, dims=state_info['qutip_dims'])
            else:
                rho = qt.Qobj(rho_data)
            
            system.set_density_matrix(rho)
        
        # 恢复哈密顿量（如果存在）
        if 'hamiltonian' in state_info:
            h_data = state_info['hamiltonian']
            
            # 如果有QuTiP dims信息，使用它创建Qobj
            if 'hamiltonian_qutip_dims' in state_info:
                h = qt.Qobj(h_data, dims=state_info['hamiltonian_qutip_dims'])
            else:
                h = qt.Qobj(h_data)
            
            system.set_hamiltonian(h)
        
        # 恢复Lindblad算符（如果存在）
        if 'lindblad_ops' in state_info and state_info['lindblad_ops']:
            lindblad_ops = []
            
            for i, op_data in enumerate(state_info['lindblad_ops']):
                # 如果有QuTiP dims信息，使用它创建Qobj
                if 'lindblad_qutip_dims' in state_info and i < len(state_info['lindblad_qutip_dims']):
                    op = qt.Qobj(op_data, dims=state_info['lindblad_qutip_dims'][i])
                else:
                    op = qt.Qobj(op_data)
                
                lindblad_ops.append(op)
            
            system.set_lindblad_ops(lindblad_ops)
        
        # 如果有任何参数信息，也可以恢复
        if 'parameters' in state_info and isinstance(state_info['parameters'], dict):
            for key, value in state_info['parameters'].items():
                system.register_parameter(key, value)
        
        print(f"已从帧 {frame_index} (t={frame['t']:.6f}) 恢复量子系统状态")
