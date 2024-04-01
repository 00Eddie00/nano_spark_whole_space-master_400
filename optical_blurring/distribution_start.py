import matplotlib.pyplot as plt
import matplotlib

from optical_blurring.distribution_generator_whole import optical_blurring
from species_blurring.cal_properties_v2 import *


def common_temporal_distribution(version, dirname, my_position_list, is_continue):
    # 是否卷积
    is_conv = True
    # 光学模糊
    optical_blurring(dirname, my_position_list, version, is_conv, is_continue)
    # 画图
    draw(version, dirname, my_position_list, is_conv)
    # 是否卷积
    is_conv = False
    # 光学模糊
    optical_blurring(dirname, my_position_list, version, is_conv, is_continue)
    # 画图
    draw(version, dirname, my_position_list, is_conv)


def draw(version, species, my_position_list, is_conv):
    dir_path = "../result/"
    # 时间间隔
    time_interval = 2 * 10 ** -6 * 100 * 1000
    # 根据是否卷积给文件名不同名称
    mark = "no_conv"
    if is_conv:
        mark = "psf"

    # 卷积画图
    dis_dirname = f"{version}"
    for i in range(len(my_position_list)):
        # 文件路径，即optical_blurring生成的文件
        dis_suffix = f"_{mark}_{version}_({my_position_list[i][0]},{my_position_list[i][1]}).csv"
        temporal_path = f"{dir_path}{dis_dirname}/{species}{dis_suffix}"
        # 画图
        temporal_plotter_no_FDHM(temporal_path, time_interval, save=True)



def pre_conv_parameter(version, my_position_list):
    # 是否从原有文件继续生成
    is_continue = False
    # 参数版本
    dirname = "CaF"
    common_temporal_distribution(version, dirname, my_position_list, is_continue)
    dirname = "CaG"
    # 生成连续空间的浓度矩阵
    common_temporal_distribution(version, dirname, my_position_list, is_continue)
