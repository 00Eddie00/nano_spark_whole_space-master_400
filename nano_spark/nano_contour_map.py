import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from optical_blurring.concentration_matrix_generator_whole import cal_nano_concentration, interpolation_calculation, \
    process_concentration
from tool.cal_bcnl import cal_elements
from tool.point_scatter_util import point_scatter
from tool.tool_mkdir import *


def gen_contour(x_coords, y_coords, values, t, d, dirname):
    # 定义绘图范围
    xmin, xmax = min(x_coords), max(x_coords)
    ymin, ymax = min(y_coords), max(y_coords)
    # 创建一个新的坐标网格
    grid_x, grid_y = np.mgrid[xmin:xmax:1000j, ymin:ymax:1000j]
    # 使用原始数据进行插值
    grid_z = griddata((x_coords, y_coords), values, (grid_x, grid_y), method='cubic')
    # 绘制等高线图
    plt.figure(figsize=(8, 6))
    contour = plt.contourf(grid_x, grid_y, grid_z, cmap='viridis')
    plt.colorbar(contour, label='Concentration')
    # # 画点
    # plt.scatter(x_coords, y_coords, c=values, cmap='viridis', edgecolors='k')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'{d} {t}ms Contour Map')
    plt.grid(True)
    plt.savefig(f"{dirname}/{d}_{t}ms.png")


def gen_line(nano_original_concentration, relations, a_arr, b_arr, c_arr, nods, t, d, dirname, grids_rz,
             open_original_concentration, coordinates_dict, xy_index, radius_list):
    concentration = np.empty(501)
    for i in range(501):
        if i <= 200:
            concentration[i] = cal_nano_concentration(i, 0, nano_original_concentration, relations, a_arr, b_arr, c_arr,
                                                      nods)
        else:
            j = 0
            k = 0
            radius = i
            height = k
            # 取出将其包含在内的网格的rz
            lower_r, upper_r, lower_z, upper_z = grids_rz[i, j, k]
            # 插值计算开放空间待插入点浓度
            concentration[i] = interpolation_calculation(lower_r, upper_r, lower_z, upper_z,
                                                         radius,
                                                         height,
                                                         open_original_concentration,
                                                         coordinates_dict, nano_original_concentration,
                                                         xy_index,
                                                         radius_list, relations, a_arr, b_arr, c_arr,
                                                         nods)
    plt.figure(figsize=(8, 6))
    # 使用matplotlib绘制折线图
    plt.plot(concentration)
    # 添加标题和标签
    plt.title(f'{d} {t}ms Line Chart')
    plt.xlabel('X Coordinate')
    plt.ylabel('Concentration')
    plt.grid(True)
    plt.savefig(f"{dirname}/{d}_{t}ms.png")


def draw_contour(version):
    folder = f"../result/{version}_nano_concentration_graph"
    name1 = f"{folder}/nano_contour_map"
    name2 = f"{folder}/line_chart"
    name3 = f"{folder}/open_contour_map"
    mkdir(f"{name1}")
    mkdir(f"{name2}")
    mkdir(f"{name3}")
    time_interval = 2 * 10 ** -6 * 100 * 1000

    nods = np.loadtxt("../config/nano/4RYRnod.dat", dtype=int) - 1
    # relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
    relations = np.load("../optical_blurring/relation/refined_relations.npy")
    # 开放空间网格点坐标
    open_grids = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    z_coords = open_grids[:, 0]
    r_coords = open_grids[:, 1]
    # 创建一个哈希表（字典）用于存储坐标值与索引的映射
    coordinates_dict = {(z, r): i for i, (z, r) in enumerate(open_grids)}
    # xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
    xy_index = np.load(f"../optical_blurring/xy_list/xy_index.npy")
    # radius_list：开放空间r方向上0~300的坐标
    new_open_r = np.flipud(point_scatter(195, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (200,)))
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        "../config/nano/4RYRgridt.dat", "../config/nano/4RYRnod.dat")
    grids_zr = np.load(f"../optical_blurring/grids_zr/grids_zr_(500,0).npy")

    dirnames = ["Ca", "CaF", "CaG"]
    # 500步对应1ms
    time_list = [0, 10, 20, 30, 40]

    # 用于gen_contour
    nano_grids = np.loadtxt("../config/nano/4RYRgridt.dat", dtype=np.float64)
    x_coords = nano_grids[:, 0]
    y_coords = nano_grids[:, 1]

    nano_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\NANO_{version}_parameters"
    open_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\OPEN_{version}_parameters"
    for d in dirnames:
        # 初始化各点Ca浓度
        nano_filenames = os.listdir(f"{nano_path}\\{d}")  # 纳米空间目前所有步数的浓度文件
        open_filenames = os.listdir(f"{open_path}\\{d}")  # 开放空间目前所有步数的浓度文件
        for t in time_list:
            print(f"{d}_{t}ms开始")
            step = int(t / time_interval)
            nano_file = nano_filenames[step]
            open_file = open_filenames[step]
            nano_original_concentration = np.loadtxt(f"{nano_path}\\{d}\\{nano_file}")
            open_original_concentration = np.loadtxt(f"{open_path}\\{d}\\{open_file}")
            gen_contour(x_coords, y_coords, nano_original_concentration, t, d, name1)
            gen_contour(r_coords, z_coords, open_original_concentration, t, d, name3)
            gen_line(nano_original_concentration, relations, a_arr, b_arr, c_arr, nods, t, d, name2, grids_zr,
                     open_original_concentration, coordinates_dict, xy_index, radius_list)
