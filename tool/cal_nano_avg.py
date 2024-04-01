import os

from matplotlib import pyplot as plt

from tool.cal_bcnl import *
from nano_spark.nano_parameters import *
import numpy as np

# 计算各个点平均值
from tool.tool_mkdir import mkdir


def cal_dye_avg(version):
    time_interval = 2 * 10 ** -6 * 100 * 1000
    grid_file_name = "../config/nano/4RYRgridt.dat"
    nod_file_name = "../config/nano/4RYRnod.dat"
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    dirnames = ["Ca", "CaG"]
    prefix = f"../result/avg/{version}"
    mkdir(f"{prefix}")
    nano_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\NANO_{version}_parameters"
    for dirname in dirnames:
        # 初始化各点Ca浓度
        filenames = os.listdir(f"{nano_path}\\{dirname}")  # 目前所有步数的浓度文件
        current_step = len(filenames)  # 当前生成文件数
        c_avg = np.zeros(current_step)
        i = 0
        for filename in filenames:
            print(f"{filename}开始")
            nano_c = np.loadtxt(f"{nano_path}\\{dirname}\\{filename}")
            total = np.sum(nano_c * control_area)
            total = total / 3.0
            avg = total / total_area
            c_avg[i] = avg
            i = i + 1
        t = [j * time_interval for j in range(i)]  # 时间，横坐标，单位ms
        # np.savetxt(f"{prefix}/{dirname}_avg.csv", c_avg)
        plt.figure(figsize=(8, 6))
        # 使用matplotlib绘制折线图
        plt.plot(t, c_avg, ls='-', lw=2)
        # 添加标题和标签
        plt.title(f'{dirname} avg')
        plt.xlabel('time')
        plt.ylabel('Concentration')
        plt.grid(True)
        plt.savefig(f"{prefix}/{dirname}_avg.png")


def cal_buffer_avg(version):
    time_interval = 2 * 10 ** -6 * 100 * 1000
    grid_file_name = "../config/nano/4RYRgridt.dat"
    nod_file_name = "../config/nano/4RYRnod.dat"
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        grid_file_name, nod_file_name)
    dirname = "CaB"
    prefix = f"../result/avg/{version}"
    mkdir(f"{prefix}")
    nano_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\NANO_{version}_parameters"
    # 初始化各点Ca浓度
    filenames = os.listdir(f"{nano_path}\\{dirname}")  # 目前所有步数的浓度文件
    current_step = len(filenames)  # 当前生成文件数
    c_avg = np.zeros((current_step, 4))
    i = 0
    for filename in filenames:
        print(f"{filename}开始")
        nano_c = np.loadtxt(f"{nano_path}\\{dirname}\\{filename}")
        for j in range(0, 4):
            nano_c_j = nano_c[:, j]
            total = np.sum(nano_c_j * control_area)
            total = total / 3.0
            avg = total / total_area
            c_avg[i, j] = avg
        i = i + 1
    t = [j * time_interval for j in range(i)]  # 时间，横坐标，单位ms
    # np.savetxt(f"{prefix}/{dirname}_avg.csv", c_avg)
    name = ["Calmodulin", "Troponin C", "SR", "SL"]
    for j in range(0, 4):
        plt.figure(figsize=(8, 6))
        # 使用matplotlib绘制折线图
        plt.plot(t, c_avg[:, j], ls='-', lw=2)
        # 添加标题和标签
        plt.title(f'{name[j]} avg')
        plt.xlabel('time')
        plt.ylabel('Concentration')
        plt.grid(True)
        plt.savefig(f"{prefix}/{dirname}_{name[j]}_avg.png")
