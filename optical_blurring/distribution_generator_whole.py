from nano_spark.nano_parameters import nano_grid_file_name, nod_file_name
from optical_blurring.concentration_matrix_generator_whole import process_concentration
from tool.cal_bcnl import cal_elements
from tool.point_scatter_util import point_scatter
from tool.tool_mkdir import *
import numpy as np

'''
    dirname：物质名
    kernel：三维卷积核矩阵83*83*83
    xy_c_val, z_c_val：用于初始化浓度矩阵
    position_list：三个观测点
    version：参数版本
    is_conv：是否卷积
    is_continue：是否继续
'''


def temporal_distribution(dirname, kernel, xy_c_val, z_c_val, position_list, version, is_conv, is_continue):
    print("temporal_distribution开始")
    nano_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\NANO_{version}_parameters"
    open_path = f"D:\\Projects\\SuYuTong\\DATA\\result\\OPEN_{version}_parameters"
    # 初始化各点Ca浓度
    nano_filenames = os.listdir(f"{nano_path}\\{dirname}")  # 纳米空间目前所有步数的浓度文件
    open_filenames = os.listdir(f"{open_path}\\{dirname}")  # 开放空间目前所有步数的浓度文件
    fluo_dir_list_len = len(nano_filenames)  # 当前生成文件数
    position_list_len = len(position_list)  # 需要检测的点个数
    # position_list_len：3，fluo_dir_list_len：文件个数。存储三个点不同时刻的值
    position_con = np.empty((position_list_len, fluo_dir_list_len))
    # grids_zr：该文件保存开放空间三维坐标正半轴1000*500*500大小的三维矩阵每个点所在的开放空间的网格，之后的浓度矩阵中的任意点均可在其中找到
    if position_list_len == 1:
        grids_zr = np.load(f"../optical_blurring/grids_zr/grids_zr_(300,300).npy")  # 300,300
    else:
        grids_zr = np.load(f"../optical_blurring/grids_zr/grids_zr_(500,0).npy")  # 500,0
    # xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
    xy_index = np.load(f"../optical_blurring/xy_list/xy_index.npy")
    # radius_list：开放空间r方向上0~300的坐标
    new_open_r = np.flipud(point_scatter(195, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (200,)))
    # 开放空间网格点坐标
    grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    # relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
    relations = np.load("../optical_blurring/relation/refined_relations.npy")
    # 创建一个哈希表（字典）用于存储坐标值与索引的映射
    coordinates_dict = {(z, r): i for i, (z, r) in enumerate(grid_coordinates)}
    # 计算三角形系数
    '''
    single_area：每个三角形面积数组
    control_area：每个点的控制面积数组
    near_triangle：每个点的相邻三角形编号
    index_in_triangle：每个点在当前三角形中的编号
    nix_multiply_l，niy_multiply_l：保存公式中nix*l，niy*l数组
    a_arr，b_arr，c_arr：保存公式中的ai bi ci
    total_area：纳米空间总面积
    '''
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        nano_grid_file_name, nod_file_name)
    # 保存纳米空间每个三角形三个点的编号
    nods = np.loadtxt(nod_file_name, dtype=int) - 1
    # 从0步开始计算
    current_step = 0
    # 若是继续运行则加载原文件内容到新数组中
    if is_continue:
        # 根据是否卷积加载不同文件
        mark = "no_conv"
        if is_conv:
            mark = "psf"
        for i in range(position_list_len):
            position_con[i] = np.loadtxt(
                f"../result/NANO_{version}_parameters/{dirname}_{mark}_{version}_({position_list[i][0]},{position_list[i][1]}).csv ")
        current_step = len(position_con[0])

    # 从当前步数到已有文件数
    for fluo_file_index in range(current_step, fluo_dir_list_len):
        nano_file = nano_filenames[fluo_file_index]
        open_file = open_filenames[fluo_file_index]
        print(f"{nano_file}开始")
        # 纳米空间原始浓度值和开放空间原始浓度值
        nano_original_concentration = np.loadtxt(f"{nano_path}\\{dirname}\\{nano_file}")
        open_original_concentration = np.loadtxt(f"{open_path}\\{dirname}\\{open_file}")
        # 遍历三个观察点
        for position_index in range(position_list_len):
            position_list_i = position_list[position_index]  # (0,0),(300,0),(500,0)
            # 返回浓度矩阵与中心点位置
            '''
            nano_original_concentration,open_original_concentration：纳米空间原始浓度值和开放空间原始浓度值
            xy_c_val：用于初始化浓度矩阵
            position_list_i：分别取0，300，500
            grids_zr：该文件保存开放空间三维坐标正半轴1000*500*500大小的三维矩阵每个点所在的开放空间的网格，之后的浓度矩阵中的任意点均可在其中找到
            xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
            radius_list：开放空间r方向上0~300的坐标
            coordinates_dict：（字典）用于存储坐标值与索引的映射
            relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
            a_arr, b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
            nods：保存纳米空间每个三角形三个点的编号
            '''
            processed_con_matrix, mid = process_concentration(nano_original_concentration,
                                                              open_original_concentration,
                                                              xy_c_val,
                                                              position_list_i, grids_zr, xy_index, radius_list,
                                                              coordinates_dict, relations, a_arr, b_arr, c_arr,
                                                              nods)
            # 取出浓度矩阵的中心点即所需的观察点
            position_con[position_index, fluo_file_index] = processed_con_matrix[mid - 1, mid - 1, mid - 1]
            # 若卷积
            if is_conv:
                # 浓度矩阵与卷积核矩阵对应点相乘后求和，即为观察点卷积值
                center_result = np.sum(processed_con_matrix * kernel)
                position_con[position_index, fluo_file_index] = center_result
    return position_con


def optical_blurring(dirname, position_list, version, is_conv, is_continue):
    # 卷积核矩阵
    kernel = np.load("../optical_blurring/kernel/kernel_3D.npy", allow_pickle=True)
    # 开放空间中没有
    if dirname == "CaG":
        C_VAL = 0.0  # 用于纳米钙火花，为第一步的值
    elif dirname == "CaF":
        C_VAL = 4.081632653061224858e-03  # 用于钙火花，为第一步的值
    # 用于初始化浓度矩阵
    xy_c_val = C_VAL
    z_c_val = C_VAL
    result_set = f"../result/{version}"
    mkdir(f"{result_set}")
    # 得到不同时刻三个点的结果
    '''
    dirname：物质名
    kernel：三维卷积核矩阵83*83*83
    xy_c_val, z_c_val：用于初始化浓度矩阵
    position_list：三个观测点
    version：参数版本
    is_conv：是否卷积
    is_continue：是否继续
    '''
    position_con = temporal_distribution(dirname, kernel, xy_c_val, z_c_val, position_list, version, is_conv,
                                         is_continue)
    # 根据是否卷积给文件名不同名称
    mark = "no_conv"
    if is_conv:
        mark = "psf"
    # 三个点保存到三个文件
    for position_index in range(len(position_con)):
        np.savetxt(
            f"{result_set}/{dirname}_{mark}_{version}_({position_list[position_index][0]},{position_list[position_index][1]}).csv "
            , position_con[position_index])
