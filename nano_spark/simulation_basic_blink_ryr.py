import matplotlib
from matplotlib import pyplot as plt

from cal_nano_element_blink_ryr import *
from nano_spark.nano_contour_map import draw_contour
from nano_spark.open.cal_open_elements_blink_ryr import *
from nano_parameters import *
from optical_blurring.distribution_start import pre_conv_parameter
from tool.cal_nano_avg import cal_dye_avg, cal_buffer_avg
from tool.tool_mkdir import *
import os
import re


def nano_spark(is_continue, total_steps, version, do_save):
    # ryr通道开放时间对应的步数
    release_step = int(RELEASE_TIME / DT)
    # 开放空间网格坐标
    open_grid_coordinates = np.loadtxt("../config/open/open_grid_coordinates.csv", delimiter=",")
    # 出流边界长度
    out_boundary_length = cal_out_boundary_length(grid)
    # 开放空间点个数
    point_count = len(open_grid_coordinates)
    # 纳米入流处仍然固定
    k_ryr = K_RYR
    # ca_jsr = CA_JSR

    # 文件路径前缀
    dir_name = "D:\\Projects\\SuYuTong\\DATA\\result"
    nano_c_ca_prefix = f"{dir_name}\\NANO_{version}_parameters\\Ca\\Ca"
    nano_c_cag_prefix = f"{dir_name}\\NANO_{version}_parameters\\CaG\\CaG"
    nano_c_caf_prefix = f"{dir_name}\\NANO_{version}_parameters\\CaF\\CaF"
    nano_c_cab_prefix = f"{dir_name}\\NANO_{version}_parameters\\CaB\\CaB"

    open_c_ca_prefix = f"{dir_name}\\OPEN_{version}_parameters\\Ca\\Ca"
    open_c_caf_prefix = f"{dir_name}\\OPEN_{version}_parameters\\CaF\\CaF"
    open_c_cab_prefix = f"{dir_name}\\OPEN_{version}_parameters\\CaB\\CaB"
    open_c_cag_prefix = f"{dir_name}\\OPEN_{version}_parameters\\CaG\\CaG"  # 虽然开放空间中没有cag，但是为了方便卷积，将其值全部设置为0.0
    '''
    nano_grid_file_name：纳米空间网格坐标文件名
    nod_file_name：纳米空间三角形编号文件名
    '''
    nano_grid_file_name, nod_file_name = "../config/nano/4RYRgridt.dat", "../config/nano/4RYRnod.dat"
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
    # 组成数组，当作一个参数，方便传参
    bcnl_elements = np.array(
        [single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr,
         c_arr, nmax, total_area], dtype=object)

    # 保存开放空间每个点的邻点
    neighbors = np.loadtxt("../config/open/open_neighbor.csv", int, delimiter=",")
    # 学长将系数*2
    # coefficients = np.load("../config/open/coefficient.npy") * 2
    coefficients = np.load("../config/open/coefficient.npy")
    # 使用blink计算的jsr中ca的浓度
    fn_jsr = np.loadtxt("../config/nano/fn_jsr.csv")

    # 从0开始
    current_step = 1

    # ***************************************************************************************
    # 开放空间中无该染料，但为了之后方便生成浓度矩阵，将其全部设为0
    open_cag = np.full(point_count, 0.0)
    # 继续运行***************************************************************************************
    if is_continue:
        path = "D:\\Projects\\SuYuTong\\DATA\\result\\"
        dir_nano = f"NANO_{version}_parameters\\"
        dir_open = f"OPEN_{version}_parameters\\"
        dirnames = ["Ca", "CaB", "CaF", "CaG"]
        # 初始化各点Ca浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[0]}\\")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        matches = re.findall(r'[0-9]', last_filename)
        # 将匹配的数字转换为整数或字符串，已知步数+1=当前步数
        current_step = int(''.join(matches)) + 1  # 如果需要提取的数字作为整数
        # 关闭RYR
        if current_step >= release_step:
            k_ryr = 0
        # 纳米空间与开放空间各点Ca浓度
        nano_f = np.loadtxt(f"{path}{dir_nano}{dirnames[0]}\\{last_filename}")
        open_f = np.loadtxt(f"{path}{dir_open}{dirnames[0]}\\{last_filename}")

        # 纳米空间与开放空间各点缓冲物浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[1]}\\")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        nano_cab = np.loadtxt(f"{path}{dir_nano}{dirnames[1]}\\{last_filename}")
        open_cab = np.loadtxt(f"{path}{dir_open}{dirnames[1]}\\{last_filename}")
        # 0123分别对应 Calmodulin、Troponin C、SR membrane、SL membrane
        nano_cab1 = nano_cab[:, 0]
        open_cab1 = open_cab[:, 0]
        nano_cab2 = nano_cab[:, 1]
        open_cab2 = open_cab[:, 1]
        nano_cab3 = nano_cab[:, 2]
        open_cab3 = open_cab[:, 2]
        nano_cab4 = nano_cab[:, 3]
        open_cab4 = open_cab[:, 3]

        # 初始化各点CaF浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[2]}\\")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        # 纳米空间与开放空间各点Caf浓度
        nano_caf = np.loadtxt(f"{path}{dir_nano}{dirnames[2]}\\{last_filename}")
        open_caf = np.loadtxt(f"{path}{dir_open}{dirnames[2]}\\{last_filename}")

        # 初始化纳米空间CaG浓度
        filenames = os.listdir(f"{path}{dir_nano}{dirnames[3]}\\")  # 目前所有步数的浓度文件
        last_filename = filenames[-1]
        nano_cag = np.loadtxt(f"{path}{dir_nano}{dirnames[3]}\\{last_filename}")
        nano_cag = set_cag_0(nano_cag)
        # 开放空间中无该染料

    # 从1开始****************************************************************************************
    else:
        # 初始化各点Ca浓度
        nano_f = np.full(NP, INITIAL_C_CA)
        open_f = np.full(point_count, INITIAL_C_CA)

        # 初始化各点CaF浓度
        INITIAL_F3 = K_F3_PLUS * INITIAL_C_CA * F3_T / (K_F3_PLUS * INITIAL_C_CA + K_F3_MINUS)
        nano_caf = np.full(NP, INITIAL_F3)
        open_caf = np.full(point_count, INITIAL_F3)

        # 初始化各点CaG浓度
        INITIAL_G = K_GCaMP6f_PLUS * INITIAL_C_CA * GCaMP6f_T / (K_GCaMP6f_PLUS * INITIAL_C_CA + K_GCaMP6f_MINUS)
        nano_cag = np.full(NP, INITIAL_G)
        nano_cag = set_cag_0(nano_cag)
        # 开放空间中无该染料

        # 初始化各点Calmodulin浓度
        INITIAL_H1 = K_Calmodulin_PLUS * INITIAL_C_CA * Calmodulin_T / (
                K_Calmodulin_PLUS * INITIAL_C_CA + K_Calmodulin_MINUS)
        nano_cab1 = np.full(NP, INITIAL_H1)
        open_cab1 = np.full(point_count, INITIAL_H1)

        # 初始化各点Troponin C浓度
        INITIAL_H2 = K_TroponinC_PLUS * INITIAL_C_CA * TroponinC_T / (
                K_TroponinC_PLUS * INITIAL_C_CA + K_TroponinC_MINUS)
        nano_cab2 = np.full(NP, INITIAL_H2)
        open_cab2 = np.full(point_count, INITIAL_H2)

        # 初始化各点SR浓度
        INITIAL_H3 = K_SR_PLUS * INITIAL_C_CA * SR_T / (K_SR_PLUS * INITIAL_C_CA + K_SR_MINUS)
        nano_cab3 = np.full(NP, INITIAL_H3)
        open_cab3 = np.full(point_count, INITIAL_H3)

        # 初始化各点SL浓度
        INITIAL_H4_NANO = K_SL_PLUS_NANO * INITIAL_C_CA * SL_T_NANO / (K_SL_PLUS_NANO * INITIAL_C_CA + K_SL_MINUS_NANO)
        nano_cab4 = np.full(NP, INITIAL_H4_NANO)
        INITIAL_H4_OPEN = K_SL_PLUS_OPEN * INITIAL_C_CA * SL_T_OPEN / (K_SL_PLUS_OPEN * INITIAL_C_CA + K_SL_MINUS_OPEN)
        open_cab4 = np.full(point_count, INITIAL_H4_OPEN)

        # 生成保存纳米空间的文件***********************************************************************************
        mkdir(f"{dir_name}/NANO_{version}_parameters/Ca")
        mkdir(f"{dir_name}/NANO_{version}_parameters/CaG")
        mkdir(f"{dir_name}/NANO_{version}_parameters/CaF")
        mkdir(f"{dir_name}/NANO_{version}_parameters/CaB")
        mkdir(f"{dir_name}/NANO_{version}_parameters/nano_avg")

        # 生成保存开放空间的文件
        mkdir(f"{dir_name}/OPEN_{version}_parameters/Ca")
        mkdir(f"{dir_name}/OPEN_{version}_parameters/CaG")
        mkdir(f"{dir_name}/OPEN_{version}_parameters/CaF")
        mkdir(f"{dir_name}/OPEN_{version}_parameters/CaB")

        # 纳米空间**********************************************************************
        # 记录初始时刻纳米空间钙离子浓度值
        ca_file_name = f"{nano_c_ca_prefix}00000000.csv"
        np.savetxt(ca_file_name, nano_f)
        print(ca_file_name, "SAVED")

        # 记录初始时刻纳米空间荧光（GCaMP6f）钙离子浓度值
        cag_file_name = f"{nano_c_cag_prefix}00000000.csv"
        np.savetxt(cag_file_name, nano_cag)
        print(cag_file_name, "SAVED")

        # 记录初始时刻纳米空间荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{nano_c_caf_prefix}00000000.csv"
        np.savetxt(caf_file_name, nano_caf)
        print(caf_file_name, "SAVED")

        # 记录初始时刻纳米空间缓冲物浓度值，纳米空间只有SR membrane、SL membrane
        cab_file_name = f"{nano_c_cab_prefix}00000000.csv"
        np.savetxt(cab_file_name, np.stack((nano_cab1, nano_cab2, nano_cab3, nano_cab4), 1))
        print(cab_file_name, "SAVED")

        # 开放空间*****************************************************************************************
        # 记录初始时刻开放空间钙离子浓度值
        ca_file_name = f"{open_c_ca_prefix}00000000.csv"
        np.savetxt(ca_file_name, open_f)
        print(ca_file_name, "SAVED")

        # 记录初始时刻开放空间荧光（GCaMP6f）钙离子浓度值，全部都是0.0
        cag_file_name = f"{open_c_cag_prefix}00000000.csv"
        np.savetxt(cag_file_name, open_cag)
        print(cag_file_name, "SAVED")

        # 记录初始时刻开放空间荧光（Fluo-3）钙离子浓度值
        caf_file_name = f"{open_c_caf_prefix}00000000.csv"
        np.savetxt(caf_file_name, open_caf)
        print(caf_file_name, "SAVED")

        # 记录初始时刻开放空间缓冲物浓度值，开放空间只有Calmodulin、Troponin C
        cab_file_name = f"{open_c_cab_prefix}00000000.csv"
        np.savetxt(cab_file_name, np.stack((open_cab1, open_cab2, open_cab3, open_cab4), 1))
        print(cab_file_name, "SAVED")
    # 开始计算***************************************************************************************************************
    for i in range(current_step, total_steps + 1):
        # 计算纳米空间*******************************************************************************************************
        # c_ca_out、c_caf_out直接取开放空间z=0.0,r=205.0处的值
        ca_jsr = fn_jsr[i]
        c_ca_out = open_f[69]
        c_caf_out = open_caf[69]
        new_nano_f, new_nano_cag, new_nano_cab1, new_nano_cab2, new_nano_cab3, new_nano_cab4 = nano_calculation_f(k_ryr,
                                                                                                                  nano_f,
                                                                                                                  nano_caf,
                                                                                                                  nano_cag,
                                                                                                                  nano_cab1,
                                                                                                                  nano_cab2,
                                                                                                                  nano_cab3,
                                                                                                                  nano_cab4,
                                                                                                                  ca_jsr,
                                                                                                                  c_ca_out,
                                                                                                                  bcnl_elements)
        '''
        k_ryr：ryr通道处钙离子释放的扩散系数
        nano_f, nano_caf, nano_cag, nano_cab3, nano_cab4：当前步数纳米空间ca、caf、cag、SR membrane、SL membrane浓度
        ca_jsr：jSR中的钙离子浓度
        c_ca_out：边界浓度
        bcnl_elements：和三角形有关的参数
        '''
        '''
        nano_f, nano_caf：当前步数纳米空间ca、caf浓度
        c_caf_out：边界浓度
        bcnl_elements：和三角形有关的参数
        '''
        new_nano_caf = nano_calculation_g2(nano_f, nano_caf, c_caf_out, bcnl_elements)

        # n + 1 步的浓度，计算出流边界平均值
        out_boundary_c_ca = np.sum(new_nano_f[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        out_boundary_c_caf = np.sum(new_nano_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        new_nano_f[84:205] = out_boundary_c_ca
        new_nano_caf[84:205] = out_boundary_c_caf

        # 计算开放空间********************************************************************************************************
        # n 步的浓度，计算出流边界平均值
        out_boundary_c_ca = np.sum(nano_f[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        out_boundary_c_caf = np.sum(nano_caf[84:205] * out_boundary_length) / np.sum(out_boundary_length)
        # 开放空间不计算交界处，直接取纳米空间计算结果
        open_f[:3] = out_boundary_c_ca
        open_caf[:3] = out_boundary_c_caf
        new_open_f, new_open_cab1, new_open_cab2, new_open_cab3, new_open_cab4 = open_calculation_f(open_f, open_caf,
                                                                                                    open_cab1,
                                                                                                    open_cab2,
                                                                                                    open_cab3,
                                                                                                    open_cab4,
                                                                                                    open_grid_coordinates,
                                                                                                    neighbors,
                                                                                                    coefficients)
        '''
         open_f, open_caf, open_cab1, open_cab2：当前步数纳米空间ca、caf浓度
         open_grid_coordinates：开放空间网格坐标
         neighbors：开放空间每个点的邻点
         coefficients：开放空间每个点的系数
        '''
        new_open_caf = open_calculation_caf(open_f, open_caf, open_grid_coordinates, neighbors, coefficients)

        # 保存文件********************************************************************************************************

        # 纳米空间******************************************************************************************
        nano_f = np.copy(new_nano_f)
        nano_cag = np.copy(new_nano_cag)
        nano_caf = np.copy(new_nano_caf)
        nano_cab1 = np.copy(new_nano_cab1)
        nano_cab2 = np.copy(new_nano_cab2)
        nano_cab3 = np.copy(new_nano_cab3)
        nano_cab4 = np.copy(new_nano_cab4)

        open_f = np.copy(new_open_f)
        new_open_cag = open_cag
        open_caf = np.copy(new_open_caf)
        open_cab1 = np.copy(new_open_cab1)
        open_cab2 = np.copy(new_open_cab2)
        open_cab3 = np.copy(new_open_cab3)
        open_cab4 = np.copy(new_open_cab4)

        if i % do_save == 0:
            # 生成编号
            length = len(str(i))
            rest_name = "0" * (8 - length) + str(i)
            # 记录该时刻钙离子浓度值
            ca_file_name = f"{nano_c_ca_prefix}{rest_name}.csv"
            np.savetxt(ca_file_name, new_nano_f)
            print(ca_file_name, "SAVED")
            # 记录该时刻荧光（GCaMP6f）钙离子浓度值
            cag_file_name = f"{nano_c_cag_prefix}{rest_name}.csv"
            np.savetxt(cag_file_name, new_nano_cag)
            print(cag_file_name, "SAVED")
            # 记录该时刻荧光（Fluo-3）钙离子浓度值
            caf_file_name = f"{nano_c_caf_prefix}{rest_name}.csv"
            np.savetxt(caf_file_name, new_nano_caf)
            print(caf_file_name, "SAVED")
            # 记录该时刻缓冲物浓度值
            cab_file_name = f"{nano_c_cab_prefix}{rest_name}.csv"
            np.savetxt(cab_file_name, np.stack((new_nano_cab1, new_nano_cab2, new_nano_cab3, new_nano_cab4), 1))
            print(cab_file_name, "SAVED")

            # 开放空间******************************************************************************************************
            # 记录初始时刻开放空间钙离子浓度值
            open_ca_file_name = f"{open_c_ca_prefix}{rest_name}.csv"
            np.savetxt(open_ca_file_name, new_open_f)
            print(open_ca_file_name, "SAVED")

            # 记录初始时刻开放空间荧光（GCaMP6f）钙离子浓度值
            open_cag_file_name = f"{open_c_cag_prefix}{rest_name}.csv"
            np.savetxt(open_cag_file_name, new_open_cag)
            print(open_cag_file_name, "SAVED")

            # 记录初始时刻开放空间荧光（Fluo-3）钙离子浓度值
            open_caf_file_name = f"{open_c_caf_prefix}{rest_name}.csv"
            np.savetxt(open_caf_file_name, new_open_caf)
            print(open_caf_file_name, "SAVED")

            # 记录初始时刻开放空间缓冲物浓度值
            open_cab_file_name = f"{open_c_cab_prefix}{rest_name}.csv"
            np.savetxt(open_cab_file_name, np.stack((new_open_cab1, new_open_cab2, new_open_cab3, new_open_cab4), 1))
            print(open_cab_file_name, "SAVED")

        # 控制ryr通道关闭*************************************************************************
        if i == release_step:
            k_ryr = 0


if __name__ == '__main__':
    matplotlib.use('TkAgg')
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    # 是否根据已有浓度继续运行
    is_continue = False
    # 执行步数
    total_steps = 20000
    do_save = 100
    # 参数版本
    version = "basic_200_(blink's_ryr_fn_K+%5)"
    nano_spark(is_continue, total_steps, version, do_save)
    draw_contour(version)
    # 观察点
    position_list = [[100, 0], [200, 0], [300, 0], [500, 0]]
    pre_conv_parameter(version, position_list)
    position_list = [[300, 300]]
    pre_conv_parameter(version, position_list)
    cal_dye_avg(version)
    cal_buffer_avg(version)
