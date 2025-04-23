# coding=utf-8
# metrics.py

"""
性能指标计算模块

该模块提供了计算和分析DLCZ量子接口性能指标的工具。它包括量子存储器性能测量、
纠缠质量评估、纠缠生成率计算以及实用的统计方法，用于对模拟结果进行评估和对比。

主要功能包括：
1. 量子存储器性能指标 - 保真度、相干时间、读取效率等
2. 光子源特性分析 - 单光子纯度、波包模式匹配等 
3. 纠缠质量测量 - 保真度、纠缠熵、Bell参数等
4. 系统效率指标 - 纠缠生成率、成功概率等
5. 统计分析工具 - 参数敏感性分析、结果可视化等
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any, Optional, Union
from scipy.optimize import curve_fit
import json
import os
import pandas as pd
from scipy.stats import sem

class PerformanceMetrics:
    """
    性能指标计算类
    
    提供了一系列方法来计算、分析和可视化DLCZ协议的性能指标。
    """
    
    def __init__(self, config=None):
        """
        初始化性能指标计算器
        
        参数:
            config: 配置对象，包含性能评估参数
        """
        self.config = config if config else {}
        
        # 指标参数
        self.target_fidelity = self.config.get('target_fidelity', 0.9)
        self.reference_values = self.config.get('reference_values', {})
        
        # 结果存储
        self.metrics_results = {}
        self.metric_history = {}
    
    def calculate_storage_metrics(self, storage_results):
        """
        计算量子存储器性能指标
        
        参数:
            storage_results: 存储过程结果字典
            
        返回:
            dict: 存储器性能指标
        """
        # 提取相关数据
        metrics = {}
        
        if 'times' in storage_results and 'coherence' in storage_results:
            times = storage_results['times']
            coherence = storage_results['coherence']
            
            # 计算相干时间 T2
            try:
                # 拟合指数衰减模型 y = A*exp(-t/T2)
                def exp_decay(t, a, t2):
                    return a * np.exp(-t / t2)
                
                # 筛选非零数据点
                mask = coherence > 1e-10
                if np.sum(mask) > 5:  # 至少需要5个有效数据点
                    popt, _ = curve_fit(exp_decay, times[mask], coherence[mask], 
                                       p0=[1.0, np.mean(times)])
                    t2_time = popt[1]
                    metrics['coherence_time'] = t2_time
                    metrics['coherence_time_us'] = t2_time * 1e6  # 微秒
            except Exception as e:
                print(f"拟合相干时间出错: {e}")
                # 使用替代方法：时间到衰减到初始值1/e的时间点
                try:
                    threshold = coherence[0] / np.e
                    idx = np.argmax(coherence < threshold)
                    if idx > 0:
                        t2_time = times[idx]
                        metrics['coherence_time'] = t2_time
                        metrics['coherence_time_us'] = t2_time * 1e6  # 微秒
                except:
                    pass
        
        # 存储保真度
        if 'storage_fidelity' in storage_results:
            metrics['storage_fidelity'] = storage_results['storage_fidelity']
            
        # 存储效率（存储态布居数）
        if 'p_storage' in storage_results:
            p_storage = storage_results['p_storage']
            metrics['storage_efficiency'] = p_storage[-1] if len(p_storage) > 0 else 0
            
        # 存储时间
        if 'times' in storage_results:
            times = storage_results['times']
            metrics['storage_duration'] = times[-1] if len(times) > 0 else 0
            metrics['storage_duration_us'] = metrics['storage_duration'] * 1e6  # 微秒
        
        return metrics
    
    def calculate_photon_metrics(self, write_results, read_results):
        """
        计算光子特性指标
        
        参数:
            write_results: 写入过程结果字典
            read_results: 读取过程结果字典
            
        返回:
            dict: 光子特性指标
        """
        metrics = {}
        
        # 信号光子发射概率
        if write_results and 'signal_photons' in write_results:
            signal_photons = write_results['signal_photons']
            metrics['signal_photon_prob'] = signal_photons[-1] if len(signal_photons) > 0 else 0
        
        # 闲置光子发射概率
        if read_results and 'idler_photons' in read_results:
            idler_photons = read_results['idler_photons']
            metrics['idler_photon_prob'] = idler_photons[-1] if len(idler_photons) > 0 else 0
        
        # 单光子纯度 g⁽²⁾(0)，通常由实验测量
        # 这里使用模拟值，理想单光子源 g⁽²⁾(0) = 0
        metrics['g2_zero'] = 0.05  # 示例值
        
        # 光子波包时间宽度
        if write_results and 'times' in write_results and 'signal_photons' in write_results:
            try:
                times = write_results['times']
                signal_rate = np.gradient(write_results['signal_photons'], times)
                
                # 找到发射率最大的点
                peak_idx = np.argmax(signal_rate)
                if peak_idx > 0:
                    # 计算半高全宽 (FWHM)
                    half_max = signal_rate[peak_idx] / 2
                    
                    # 寻找左半高点
                    left_idx = 0
                    for i in range(peak_idx, 0, -1):
                        if signal_rate[i] < half_max:
                            left_idx = i
                            break
                    
                    # 寻找右半高点
                    right_idx = len(times) - 1
                    for i in range(peak_idx, len(times)):
                        if signal_rate[i] < half_max:
                            right_idx = i
                            break
                    
                    # 计算FWHM
                    if left_idx < right_idx:
                        fwhm = times[right_idx] - times[left_idx]
                        metrics['photon_pulse_width'] = fwhm
                        metrics['photon_pulse_width_ns'] = fwhm * 1e9  # 纳秒
            except Exception as e:
                print(f"计算光子脉冲宽度出错: {e}")
        
        return metrics
    
    def calculate_entanglement_metrics(self, verification_results):
        """
        计算纠缠质量指标
        
        参数:
            verification_results: 纠缠验证结果字典
            
        返回:
            dict: 纠缠质量指标
        """
        metrics = {}
        
        # 纠缠保真度
        if verification_results and 'entanglement_fidelity' in verification_results:
            metrics['entanglement_fidelity'] = verification_results['entanglement_fidelity']
        
        # Bell参数
        if verification_results and 'bell_parameter' in verification_results:
            bell_parameter = verification_results['bell_parameter']
            metrics['bell_parameter'] = bell_parameter
            
            # 计算Bell不等式违背程度
            # 经典极限为2，量子力学预测的最大值为2√2≈2.83
            classical_limit = 2.0
            theoretical_max = 2.0 * np.sqrt(2)
            
            if bell_parameter > classical_limit:
                # 量子化程度 = (S - 2) / (2√2 - 2)
                quantum_degree = (bell_parameter - classical_limit) / (theoretical_max - classical_limit)
                metrics['quantum_degree'] = quantum_degree
                
                # 额外的描述性标签
                if quantum_degree < 0.33:
                    metrics['entanglement_quality'] = 'weak'
                elif quantum_degree < 0.67:
                    metrics['entanglement_quality'] = 'moderate'
                else:
                    metrics['entanglement_quality'] = 'strong'
            else:
                metrics['quantum_degree'] = 0.0
                metrics['entanglement_quality'] = 'classical'
        
        # 符合率
        if verification_results and 'coincidence_rate' in verification_results:
            metrics['coincidence_rate'] = verification_results['coincidence_rate']
        
        # 成功概率
        if verification_results and 'success_probability' in verification_results:
            metrics['success_probability'] = verification_results['success_probability']
            
            # 计算纠缠生成率 = 成功概率 * 重复率
            repetition_rate = 1e6  # 假设重复率为1MHz
            metrics['entanglement_generation_rate'] = metrics['success_probability'] * repetition_rate
        
        return metrics
    
    def calculate_dlcz_protocol_metrics(self, protocol_results):
        """
        计算DLCZ协议整体性能指标
        
        参数:
            protocol_results: 协议执行结果字典
            
        返回:
            dict: 协议性能指标
        """
        # 初始化综合指标
        metrics = {
            'protocol_success': False,
            'total_duration': 0.0,
            'entanglement_fidelity': 0.0,
            'success_probability': 0.0,
            'entanglement_generation_rate': 0.0
        }
        
        # 提取协议状态
        if 'protocol_state' in protocol_results:
            state = protocol_results['protocol_state']
            metrics['protocol_success'] = state.get('entanglement_verified', False)
            metrics['total_duration'] = state.get('elapsed_time', 0.0)
            metrics['entanglement_fidelity'] = state.get('entanglement_fidelity', 0.0)
            metrics['success_probability'] = state.get('success_probability', 0.0)
        
        # 计算各阶段指标
        stage_metrics = {}
        
        # 存储指标
        if 'stage_results' in protocol_results and 'storage' in protocol_results['stage_results']:
            storage_results = protocol_results['stage_results']['storage']
            if storage_results:
                storage_metrics = self.calculate_storage_metrics(storage_results)
                stage_metrics['storage'] = storage_metrics
                
                # 将部分存储指标提升到顶层
                if 'coherence_time_us' in storage_metrics:
                    metrics['coherence_time_us'] = storage_metrics['coherence_time_us']
                if 'storage_fidelity' in storage_metrics:
                    metrics['storage_fidelity'] = storage_metrics['storage_fidelity']
        
        # 光子指标
        write_results = None
        read_results = None
        
        if 'stage_results' in protocol_results:
            if 'write' in protocol_results['stage_results']:
                write_results = protocol_results['stage_results']['write']
            if 'read' in protocol_results['stage_results']:
                read_results = protocol_results['stage_results']['read']
                
        if write_results or read_results:
            photon_metrics = self.calculate_photon_metrics(write_results, read_results)
            stage_metrics['photon'] = photon_metrics
            
            # 将部分光子指标提升到顶层
            if 'signal_photon_prob' in photon_metrics:
                metrics['signal_photon_prob'] = photon_metrics['signal_photon_prob']
            if 'idler_photon_prob' in photon_metrics:
                metrics['idler_photon_prob'] = photon_metrics['idler_photon_prob']
        
        # 纠缠指标
        if 'stage_results' in protocol_results and 'verification' in protocol_results['stage_results']:
            verification_results = protocol_results['stage_results']['verification']
            if verification_results:
                entanglement_metrics = self.calculate_entanglement_metrics(verification_results)
                stage_metrics['entanglement'] = entanglement_metrics
                
                # 将部分纠缠指标提升到顶层
                if 'bell_parameter' in entanglement_metrics:
                    metrics['bell_parameter'] = entanglement_metrics['bell_parameter']
                if 'entanglement_generation_rate' in entanglement_metrics:
                    metrics['entanglement_generation_rate'] = entanglement_metrics['entanglement_generation_rate']
        
        # 将阶段指标添加到综合指标
        metrics['stage_metrics'] = stage_metrics
        
        # 计算协议效率
        # 协议效率 = 成功概率 * 纠缠保真度 / 协议总时间
        if metrics['total_duration'] > 0:
            protocol_efficiency = (metrics['success_probability'] * 
                                  metrics['entanglement_fidelity'] / 
                                  metrics['total_duration'])
            metrics['protocol_efficiency'] = protocol_efficiency
        
        # 存储结果
        self.metrics_results = metrics
        
        return metrics
    
    def compare_with_reference(self, metrics=None):
        """
        与参考值比较
        
        参数:
            metrics: 要比较的指标，如果为None则使用最新计算的指标
            
        返回:
            dict: 比较结果
        """
        if metrics is None:
            metrics = self.metrics_results
            
        if not metrics or not self.reference_values:
            return {}
            
        comparison = {}
        
        for key, ref_value in self.reference_values.items():
            if key in metrics:
                current_value = metrics[key]
                difference = current_value - ref_value
                percent_diff = (difference / ref_value) * 100 if ref_value != 0 else float('inf')
                
                comparison[key] = {
                    'current': current_value,
                    'reference': ref_value,
                    'difference': difference,
                    'percent_diff': percent_diff,
                    'improved': difference > 0 if 'fidelity' in key or 'rate' in key or 'time' in key 
                                else difference < 0
                }
                
        return comparison
    
    def analyze_parameter_sensitivity(self, parameter_name, parameter_values, 
                                    target_metric, simulation_func):
        """
        分析参数敏感性
        
        运行多次模拟，改变某个参数的值，观察目标指标的变化
        
        参数:
            parameter_name: 参数名称
            parameter_values: 参数值列表
            target_metric: 目标指标名称
            simulation_func: 模拟函数，接受参数值并返回模拟结果
            
        返回:
            dict: 参数敏感性分析结果
        """
        results = []
        metric_values = []
        
        for value in parameter_values:
            # 运行模拟
            sim_result = simulation_func(value)
            
            # 计算指标
            metrics = self.calculate_dlcz_protocol_metrics(sim_result)
            
            # 提取目标指标
            metric_value = metrics.get(target_metric)
            
            if metric_value is not None:
                results.append({
                    'parameter_value': value,
                    'metric_value': metric_value,
                    'simulation_result': sim_result
                })
                metric_values.append(metric_value)
        
        # 计算敏感度 (metric变化量 / parameter变化量)
        sensitivities = []
        for i in range(1, len(parameter_values)):
            delta_param = parameter_values[i] - parameter_values[i-1]
            delta_metric = metric_values[i] - metric_values[i-1]
            
            if delta_param != 0:
                sensitivity = delta_metric / delta_param
                sensitivities.append(sensitivity)
        
        # 找出最敏感的区域
        if sensitivities:
            max_sens_idx = np.argmax(np.abs(sensitivities))
            max_sensitivity = sensitivities[max_sens_idx]
            max_sens_param_value = (parameter_values[max_sens_idx] + parameter_values[max_sens_idx+1]) / 2
        else:
            max_sensitivity = 0
            max_sens_param_value = 0
        
        # 构建结果
        analysis_result = {
            'parameter_name': parameter_name,
            'parameter_values': parameter_values,
            'target_metric': target_metric,
            'metric_values': metric_values,
            'sensitivities': sensitivities,
            'max_sensitivity': max_sensitivity,
            'max_sensitivity_parameter_value': max_sens_param_value,
            'detailed_results': results
        }
        
        # 保存到历史
        self.metric_history[f"{parameter_name}_{target_metric}"] = analysis_result
        
        return analysis_result
    
    def optimize_parameters(self, parameter_ranges, target_metric, simulation_func, 
                          optimization_method='grid', num_points=5):
        """
        优化参数
        
        寻找能够最大化或最小化目标指标的参数组合
        
        参数:
            parameter_ranges: 参数范围字典，格式为{参数名: (最小值, 最大值)}
            target_metric: 目标指标名称
            simulation_func: 模拟函数，接受参数字典并返回模拟结果
            optimization_method: 优化方法，'grid'或'random'
            num_points: 每个参数的采样点数（对于grid方法）或总采样点数（对于random方法）
            
        返回:
            dict: 优化结果
        """
        results = []
        
        if optimization_method == 'grid':
            # 网格搜索
            # 为每个参数生成网格点
            grid_points = {}
            for param, (min_val, max_val) in parameter_ranges.items():
                grid_points[param] = np.linspace(min_val, max_val, num_points)
            
            # 生成参数组合
            import itertools
            param_combinations = []
            param_names = list(parameter_ranges.keys())
            for values in itertools.product(*[grid_points[param] for param in param_names]):
                param_dict = {param_names[i]: values[i] for i in range(len(param_names))}
                param_combinations.append(param_dict)
        else:
            # 随机搜索
            param_combinations = []
            for _ in range(num_points):
                param_dict = {}
                for param, (min_val, max_val) in parameter_ranges.items():
                    param_dict[param] = min_val + np.random.random() * (max_val - min_val)
                param_combinations.append(param_dict)
        
        # 评估每个参数组合
        for params in param_combinations:
            # 运行模拟
            sim_result = simulation_func(params)
            
            # 计算指标
            metrics = self.calculate_dlcz_protocol_metrics(sim_result)
            
            # 提取目标指标
            metric_value = metrics.get(target_metric)
            
            if metric_value is not None:
                results.append({
                    'parameters': params,
                    'metric_value': metric_value,
                    'metrics': metrics
                })
        
        # 按目标指标排序
        # 假设我们要最大化指标（如保真度）
        results.sort(key=lambda x: x['metric_value'], reverse=True)
        
        # 构建优化结果
        optimization_result = {
            'best_parameters': results[0]['parameters'] if results else None,
            'best_metric_value': results[0]['metric_value'] if results else None,
            'all_results': results,
            'target_metric': target_metric,
            'parameter_ranges': parameter_ranges,
            'optimization_method': optimization_method
        }
        
        return optimization_result
    
    def save_metrics(self, filename, metrics=None):
        """
        保存指标到文件
        
        参数:
            filename: 文件名
            metrics: 要保存的指标，如果为None则使用最新计算的指标
            
        返回:
            bool: 保存是否成功
        """
        if metrics is None:
            metrics = self.metrics_results
            
        if not metrics:
            return False
            
        try:
            # 确保目录存在
            os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
            
            # 保存为JSON
            with open(filename, 'w') as f:
                json.dump(metrics, f, indent=2)
                
            return True
        except Exception as e:
            print(f"保存指标时出错: {e}")
            return False
    
    def load_metrics(self, filename):
        """
        从文件加载指标
        
        参数:
            filename: 文件名
            
        返回:
            dict: 加载的指标
        """
        try:
            with open(filename, 'r') as f:
                metrics = json.load(f)
                
            self.metrics_results = metrics
            return metrics
        except Exception as e:
            print(f"加载指标时出错: {e}")
            return {}
    
    def plot_metrics(self, metrics=None, figsize=(12, 8)):
        """
        绘制指标图表
        
        参数:
            metrics: 要绘制的指标，如果为None则使用最新计算的指标
            figsize: 图像大小
            
        返回:
            plt.Figure: 图像对象
        """
        if metrics is None:
            metrics = self.metrics_results
            
        if not metrics:
            return None
            
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        
        # 1. 纠缠指标
        ax = axes[0]
        entanglement_metrics = [
            metrics.get('entanglement_fidelity', 0),
            metrics.get('bell_parameter', 0) / 3,  # 归一化Bell参数
            metrics.get('quantum_degree', 0) if 'quantum_degree' in metrics else 0
        ]
        labels = ['Fidelity', 'Bell Parameter/3', 'Quantum Degree']
        
        ax.bar(labels, entanglement_metrics)
        ax.set_ylim(0, 1.1)
        ax.set_title('Entanglement Metrics')
        ax.grid(True, alpha=0.3)
        
        # 2. 存储性能
        ax = axes[1]
        coherence_time = metrics.get('coherence_time_us', 0)
        storage_fidelity = metrics.get('storage_fidelity', 0)
        
        ax.bar(['Coherence Time (μs)', 'Storage Fidelity'], [coherence_time, storage_fidelity])
        ax.set_title('Storage Performance')
        ax.grid(True, alpha=0.3)
        
        # 3. 效率指标
        ax = axes[2]
        efficiency_metrics = [
            metrics.get('success_probability', 0),
            metrics.get('signal_photon_prob', 0),
            metrics.get('idler_photon_prob', 0)
        ]
        labels = ['Success Prob.', 'Signal Photon', 'Idler Photon']
        
        ax.bar(labels, efficiency_metrics)
        ax.set_ylim(0, 1.1)
        ax.set_title('Efficiency Metrics')
        ax.grid(True, alpha=0.3)
        
        # 4. 生成率
        ax = axes[3]
        generation_rate = metrics.get('entanglement_generation_rate', 0)
        
        ax.bar(['Entanglement Generation Rate (Hz)'], [generation_rate])
        ax.set_title('Generation Rate')
        ax.grid(True, alpha=0.3)
        
        # 添加总标题
        protocol_success = metrics.get('protocol_success', False)
        fig.suptitle(f'DLCZ Protocol Performance Metrics - {"SUCCESS" if protocol_success else "FAILED"}',
                   fontsize=14, y=0.98)
        
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
    
    def plot_parameter_sensitivity(self, analysis_result, figsize=(10, 6)):
        """
        绘制参数敏感性分析图表
        
        参数:
            analysis_result: 参数敏感性分析结果
            figsize: 图像大小
            
        返回:
            plt.Figure: 图像对象
        """
        if not analysis_result:
            return None
            
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # 1. 绘制指标与参数的关系
        ax = axes[0]
        parameter_values = analysis_result['parameter_values']
        metric_values = analysis_result['metric_values']
        
        ax.plot(parameter_values, metric_values, 'bo-')
        ax.set_xlabel(analysis_result['parameter_name'])
        ax.set_ylabel(analysis_result['target_metric'])
        ax.set_title(f'{analysis_result["target_metric"]} vs {analysis_result["parameter_name"]}')
        ax.grid(True)
        
        # 2. 绘制敏感度
        ax = axes[1]
        sensitivities = analysis_result['sensitivities']
        sensitivity_x = [(parameter_values[i] + parameter_values[i+1])/2 for i in range(len(parameter_values)-1)]
        
        ax.plot(sensitivity_x, sensitivities, 'ro-')
        ax.set_xlabel(analysis_result['parameter_name'])
        ax.set_ylabel(f'Sensitivity of {analysis_result["target_metric"]}')
        ax.set_title('Parameter Sensitivity')
        ax.grid(True)
        
        # 标记最敏感点
        if sensitivities:
            max_sens_idx = np.argmax(np.abs(sensitivities))
            max_sensitivity = sensitivities[max_sens_idx]
            max_sens_param_value = sensitivity_x[max_sens_idx]
            
            ax.plot(max_sens_param_value, max_sensitivity, 'go', markersize=10)
            ax.annotate(f'Max Sensitivity: {max_sensitivity:.2e}', 
                      xy=(max_sens_param_value, max_sensitivity), 
                      xytext=(0, 20), textcoords='offset points',
                      ha='center', arrowprops=dict(arrowstyle="->"))
        
        fig.tight_layout()
        return fig
    
    def plot_optimization_results(self, optimization_result, figsize=(12, 10)):
        """
        绘制参数优化结果图表
        
        参数:
            optimization_result: 参数优化结果
            figsize: 图像大小
            
        返回:
            plt.Figure: 图像对象
        """
        if not optimization_result or not optimization_result['all_results']:
            return None
            
        # 提取数据
        results = optimization_result['all_results']
        target_metric = optimization_result['target_metric']
        best_params = optimization_result['best_parameters']
        best_metric = optimization_result['best_metric_value']
        
        # 如果有过多结果，只显示前20个
        if len(results) > 20:
            results = results[:20]
        
        # 创建DataFrame方便处理和绘图
        df = pd.DataFrame([r['parameters'] for r in results])
        df['metric_value'] = [r['metric_value'] for r in results]
        
        # 获取参数名称
        param_names = list(best_params.keys())
        
        # 绘图数量取决于参数数量
        n_params = len(param_names)
        n_rows = min(2, n_params)
        n_cols = int(np.ceil(n_params / n_rows)) + 1  # 加1是为了指标柱状图
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        axes = axes.flatten()
        
        # 绘制每个参数与指标的关系
        for i, param in enumerate(param_names):
            if i < len(axes):
                ax = axes[i]
                ax.scatter(df[param], df['metric_value'])
                ax.set_xlabel(param)
                ax.set_ylabel(target_metric)
                ax.set_title(f'{target_metric} vs {param}')
                ax.grid(True)
                
                # 标记最优点
                ax.plot(best_params[param], best_metric, 'ro', markersize=10)
        
        # 绘制最佳指标值
        ax = axes[-1]
        ax.bar(['Best ' + target_metric], [best_metric])
        ax.set_title(f'Best {target_metric}: {best_metric:.4f}')
        ax.grid(True)
        
        # 添加总标题
        param_str = ', '.join([f'{k}={v:.3g}' for k, v in best_params.items()])
        fig.suptitle(f'Parameter Optimization Results\nBest Parameters: {param_str}',
                   fontsize=14, y=0.98)
        
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
    
    def analyze_experiment_results(self, experiment_results, group_by=None, confidence=0.95):
        """
        分析实验结果
        
        参数:
            experiment_results: 实验结果列表，每个元素为一次实验的指标字典
            group_by: 分组参数名称
            confidence: 置信区间
            
        返回:
            dict: 实验分析结果
        """
        if not experiment_results:
            return {}
            
        # 提取常见指标
        common_metrics = set.intersection(*[set(res.keys()) for res in experiment_results])
        
        # 基本统计分析
        stats = {}
        for metric in common_metrics:
            values = [res[metric] for res in experiment_results if isinstance(res[metric], (int, float))]
            if values:
                stats[metric] = {
                    'mean': np.mean(values),
                    'median': np.median(values),
                    'std': np.std(values),
                    'min': np.min(values),
                    'max': np.max(values),
                    'count': len(values)
                }
                
                # 计算置信区间
                sem_val = sem(values)
                from scipy import stats as scipy_stats
                ci = scipy_stats.t.interval(
                    confidence, len(values)-1, loc=np.mean(values), scale=sem_val)
                stats[metric]['ci_lower'] = ci[0]
                stats[metric]['ci_upper'] = ci[1]
        
        # 分组分析
        grouped_stats = {}
        if group_by is not None:
            # 获取该参数的所有不同值
            group_values = list(set(res.get(group_by) for res in experiment_results 
                                 if group_by in res))
            
            for group_val in group_values:
                # 筛选该组的结果
                group_results = [res for res in experiment_results 
                               if group_by in res and res[group_by] == group_val]
                
                # 计算该组的统计
                group_stats = {}
                for metric in common_metrics:
                    values = [res[metric] for res in group_results 
                             if isinstance(res[metric], (int, float))]
                    if values:
                        group_stats[metric] = {
                            'mean': np.mean(values),
                            'std': np.std(values),
                            'count': len(values)
                        }
                
                grouped_stats[group_val] = group_stats
        
        # 构建结果
        analysis = {
            'stats': stats,
            'grouped_stats': grouped_stats,
            'sample_size': len(experiment_results),
            'common_metrics': list(common_metrics),
            'group_by': group_by
        }
        
        return analysis
    
    def plot_experiment_analysis(self, analysis, metrics_to_plot=None, figsize=(12, 10)):
        """
        绘制实验分析图表
        
        参数:
            analysis: 实验分析结果
            metrics_to_plot: 要绘制的指标列表，如果为None则使用所有指标
            figsize: 图像大小
            
        返回:
            plt.Figure: 图像对象
        """
        if not analysis or 'stats' not in analysis:
            return None
            
        # 确定要绘制的指标
        if metrics_to_plot is None:
            metrics_to_plot = list(analysis['stats'].keys())[:4]  # 默认绘制前4个指标
        
        # 计算所需子图数量
        n_metrics = len(metrics_to_plot)
        n_cols = min(2, n_metrics)
        n_rows = int(np.ceil(n_metrics / n_cols))
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        if n_rows == 1 and n_cols == 1:
            axes = np.array([axes])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.reshape(-1)
        
        for i, metric in enumerate(metrics_to_plot[:n_rows*n_cols]):
            ax = axes[i // n_cols, i % n_cols] if n_rows > 1 and n_cols > 1 else axes[i]
            
            # 获取统计数据
            stats = analysis['stats'].get(metric, {})
            if not stats:
                continue
                
            # 如果有分组数据，绘制分组比较
            if analysis['group_by'] and 'grouped_stats' in analysis:
                groups = []
                means = []
                stds = []
                
                for group_val, group_stats in analysis['grouped_stats'].items():
                    if metric in group_stats:
                        groups.append(str(group_val))
                        means.append(group_stats[metric]['mean'])
                        stds.append(group_stats[metric]['std'])
                
                if groups:
                    ax.bar(groups, means, yerr=stds, capsize=5)
                    ax.set_title(f'{metric} by {analysis["group_by"]}')
                    ax.grid(True, alpha=0.3)
                    # 旋转x轴标签如果太多
                    if len(groups) > 5:
                        ax.set_xticklabels(groups, rotation=45, ha='right')
            else:
                # 绘制单个指标的统计
                ax.bar(['Mean Value'], [stats['mean']], yerr=[stats['std']], capsize=10)
                ax.errorbar(['Mean Value'], [stats['mean']], 
                          yerr=[[stats['mean'] - stats['ci_lower']], [stats['ci_upper'] - stats['mean']]],
                          fmt='ro', capsize=5, label=f'{int(analysis.get("confidence", 0.95)*100)}% CI')
                
                ax.set_title(f'{metric} Statistics')
                ax.grid(True, alpha=0.3)
                ax.legend()
                
                # 添加关键统计
                ax.text(0.05, 0.95, f"Mean: {stats['mean']:.4f}\nStd: {stats['std']:.4f}\nN: {stats['count']}",
                      transform=ax.transAxes, va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # 添加总标题
        fig.suptitle(f'Experiment Analysis (N={analysis["sample_size"]})',
                   fontsize=14, y=0.98)
        
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
        
# 辅助函数
def exponential_decay(t, a, tau):
    """指数衰减函数: a * exp(-t/tau)"""
    return a * np.exp(-t / tau)

def calculate_fidelity_threshold(error_rate, distance, loss_alpha):
    """
    计算给定距离下的纠缠保真度阈值
    
    参数:
        error_rate: 量子比特误差率
        distance: 距离 (km)
        loss_alpha: 光纤损耗系数 (dB/km)
        
    返回:
        float: 保真度阈值
    """
    # 计算信道透射率
    transmittance = 10**(-loss_alpha * distance / 10)
    
    # 简化的保真度阈值模型
    F_threshold = 0.5 + 0.5 * transmittance * (1 - 2*error_rate)
    
    return min(1.0, max(0.5, F_threshold))

def calculate_bell_inequality(rho, method="CHSH"):
    """
    计算Bell不等式违背程度
    
    参数:
        rho: 量子态密度矩阵
        method: Bell不等式类型，支持"CHSH"(默认)
        
    返回:
        float: Bell参数，对CHSH，经典极限为2，量子极限为2*sqrt(2)≈2.82
    """
    from qutip import expect, sigmax, sigmay, tensor
    
    try:
        # 确保是密度矩阵
        if hasattr(rho, 'type') and rho.type != 'oper':
            raise ValueError("输入必须是密度矩阵（qutip.Qobj，类型为'oper'）")
        
        # 检查子系统数量
        if len(rho.dims[0]) != 2:
            raise ValueError(f"仅支持两个子系统的Bell测试，当前子系统数量: {len(rho.dims[0])}")
        
        # 检查子系统维度
        if rho.dims[0][0] != rho.dims[0][1] or rho.dims[0][0] != 2:
            raise ValueError(f"Bell测试要求两个2维子系统，当前维度: {rho.dims[0]}")
        
        # 设置测量算符
        sigma_x1 = tensor(sigmax(), np.eye(2))
        sigma_z1 = tensor(np.array([[1, 0], [0, -1]]), np.eye(2))
        sigma_x2 = tensor(np.eye(2), sigmax())
        sigma_z2 = tensor(np.eye(2), np.array([[1, 0], [0, -1]]))
        
        # 对于CHSH不等式
        # S = ⟨A₁B₁⟩ + ⟨A₁B₂⟩ + ⟨A₂B₁⟩ - ⟨A₂B₂⟩
        # 其中A₁=X₁, A₂=Z₁, B₁=(X₂+Z₂)/√2, B₂=(X₂-Z₂)/√2
        
        if method == "CHSH":
            # 测量角度
            a1 = sigma_x1
            a2 = sigma_z1
            b1 = (sigma_x2 + sigma_z2) / np.sqrt(2)
            b2 = (sigma_x2 - sigma_z2) / np.sqrt(2)
            
            # 计算关联函数
            corr11 = expect(a1 * b1, rho)
            corr12 = expect(a1 * b2, rho)
            corr21 = expect(a2 * b1, rho)
            corr22 = expect(a2 * b2, rho)
            
            # 计算CHSH参数
            S = corr11 + corr12 + corr21 - corr22
            
            return abs(S)
        else:
            raise ValueError(f"不支持的Bell不等式类型: {method}")
    
    except Exception as e:
        print(f"计算Bell不等式时出错: {e}")
        import traceback
        traceback.print_exc()
        return 0.0 