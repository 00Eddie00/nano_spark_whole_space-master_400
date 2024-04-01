import numpy as np

from tool.point_scatter_util import point_scatter

"""
arr:0~5000/0~2500
coordinate:0~1118/0~500
"""


def cal_pre_next(arr, coordinate):
    index = np.searchsorted(arr, coordinate)
    if index == 0:
        lower, upper = arr[index], arr[index + 1]
    elif index == len(arr):
        lower, upper = arr[index - 1], arr[index]
    else:
        lower, upper = arr[index - 1], arr[index]
    return lower, upper


# 只生成在三维坐标正半轴上的区域
def open_judge_relation(my_position_list):
    open_z = point_scatter(13, 2500, 6)
    open_z.insert(0, 7.5)
    open_z.insert(0, 3)
    open_z.insert(0, 0)
    open_r = point_scatter(200, 5000, 5, k=2)
    new_open_r = np.flipud(point_scatter(195, 0, 5, k=2, positive=False))
    r = np.concatenate((new_open_r, open_r))
    z = np.copy(open_z)
    # (0,0)因为只取到距离中心点400nm（0，400）处，且浓度体积为1立方微米，所以r:-500~500，z:-500~500
    # 总共要取[[0, 0], [100, 0], [300, 0], [400, 0]]
    for position_index in range(len(my_position_list)):
        dis_r = my_position_list[position_index][0]  # 300
        dis_z = my_position_list[position_index][1]  # 300
        x_end = dis_r + 500  # 800
        y_end = 500
        z_end = dis_z + 500  # 800
        # 每隔1nm取一个点
        interval = 1
        # 0~800 的平方的数组
        x_squared = np.square(np.arange(0, x_end + 1, interval))
        # 0~500 的平方的数组
        y_squared = np.square(np.arange(0, y_end + 1, interval))
        # 0~800 的数组
        z_arange = np.arange(0, z_end + 1, interval)
        d1, d2, d3 = (x_end - 0) // interval + 1, (y_end - 0) // interval + 1, (
                z_end - 0) // interval + 1  # 801 501 801
        # 待插入点的相邻点
        grids_zr = np.full((d1, d2, d3, 4), fill_value=-1.0)
        # 801  x上变化
        for i in range(d1):
            x2 = x_squared[i]
            print(f"{i}点开始")
            # 501  y上变化
            for j in range(d2):
                y2 = y_squared[j]
                radius = np.sqrt(x2 + y2)
                # 找到待插入点r的位置的前后坐标
                lower_r, upper_r = cal_pre_next(r, radius)
                # 801  z上变化
                for k in range(d3):
                    # z的坐标就是待插入位置
                    height = z_arange[k]
                    # 待插入点不在纳米空间内就插入该点
                    if height >= 7.5 or radius >= 200:
                        # 找到待插入点z的位置的前后坐标
                        lower_z, upper_z = cal_pre_next(z, height)
                        grids_zr[i, j, k] = lower_r, upper_r, lower_z, upper_z
        np.save(
            f"../optical_blurring/grids_zr/grids_zr_({my_position_list[position_index][0]},{my_position_list[position_index][1]})",
            grids_zr)


def main():
    my_position_list = [[500, 0]]
    open_judge_relation(my_position_list)
    my_position_list = [[300, 300]]
    open_judge_relation(my_position_list)
    # grids_zr = np.load("../optical_blurring/grids_zr/grids_zr_v3_(300,300).npy")
    # print(grids_zr)


if __name__ == "__main__":
    main()
