import numpy as np

grid_coordinates = np.loadtxt("open_grid_coordinates.csv", delimiter=",")
point_count = len(grid_coordinates)
big_open_count = 4692
# 建立一个二维数组存储各点的邻居；第一维表示点编号，第二维表示四个方向的此点的邻居
# 索引0：z轴正方向，索引1：z轴负方向，索引2：ρ轴正方向(外侧面)，索引3：ρ负方向(内侧面)；
# -1 代表此方向上无邻居。
neighbor = np.full((point_count, 4), -1)

# 亚空间与大开放空间交界点
boundary_r = 200
boundary_z = [0, 3, 7.5]
# 大开放区域
big_z_min = 0
big_z_max = 2500
big_r_min = 200
big_r_max = 5000
big_z_length = 69  # z轴点的个数
for i in range(0, big_open_count):
    # z轴邻居
    z = grid_coordinates[i][0]
    if big_z_min < z < big_z_max:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif big_z_min == z:
        neighbor[i][0] = i + 1
        # 以 z = 0 对称，
        # 存在下顶面，只是下顶面与上顶面浓度、面积相同，而z轴坐标互为相反数；
        # 此处将下顶面编号赋值为上顶面编号，方便取浓度值，但是计算z轴坐标相关数值时需要注意取反。
        neighbor[i][1] = i + 1
    elif big_z_max == z:
        neighbor[i][1] = i - 1
    else:
        print(f"结点 {i} z轴坐标错误")
    # r轴邻居
    r = grid_coordinates[i][1]
    if big_r_min < r < big_r_max:
        neighbor[i][2] = i + big_z_length
        neighbor[i][3] = i - big_z_length
    elif big_r_min == r:  # 与亚空间、小开放空间的交界处
        neighbor[i][2] = i + big_z_length
        # 与小开放空间交界
        if z not in boundary_z:
            neighbor[i][3] = i + 5481
    elif big_r_max == r:
        neighbor[i][3] = i - big_z_length
    else:
        print(f"结点 {i} r轴坐标错误")

# 小开放区域
small_z_min = 13
small_z_max = 2500
small_r_min = 0
small_r_max = 195
small_z_length = 66
for i in range(big_open_count, point_count):
    # z轴邻居
    z = grid_coordinates[i][0]
    if small_z_min < z < small_z_max:
        neighbor[i][0] = i + 1
        neighbor[i][1] = i - 1
    elif small_z_min == z:
        neighbor[i][0] = i + 1
    elif small_z_max == z:
        neighbor[i][1] = i - 1
    else:
        print(f"结点 {i} z轴坐标错误")
    # r轴邻居
    r = grid_coordinates[i][1]
    if small_r_min < r < small_r_max:
        neighbor[i][2] = i + small_z_length
        neighbor[i][3] = i - small_z_length
    elif small_r_min == r:
        neighbor[i][2] = i + small_z_length
    elif small_r_max == r:
        neighbor[i][2] = i - 5481
        neighbor[i][3] = i - small_z_length
    else:
        print(f"结点 {i} r轴坐标错误")

np.savetxt("open_neighbor.csv", neighbor, fmt='%d', delimiter=",")
