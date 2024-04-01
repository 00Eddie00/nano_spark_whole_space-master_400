from nano_spark.nano_parameters import *
from tool.cal_bcnl import cal_elements
from scipy.ndimage import convolve1d
from tool.generate_nano_con import *
from optical_blurring.ob_parameters import half_length


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    c_val = xy_c_val
    i = 0
    print("convolve3d开始")
    for kernel_i in kernel:
        if i == 2:
            c_val = z_c_val
        # 在第一个维度上进行一维卷积操作
        original_matrix = convolve1d(original_matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        i = i + 1
    return original_matrix


def nano_process_concentration(original_concentration, c_val, relations):
    concentration = np.full((401, 401, 16), c_val, dtype=float)
    nods = np.loadtxt(nod_file_name, dtype=int) - 1
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = cal_elements(
        nano_grid_file_name, nod_file_name)
    print("nano_process_concentration开始")
    for i in range(401):
        for j in range(401):
            x = i - 200
            y = j - 200
            relation = relations[i, j, 0]
            c_id = relations[i, j, 1]
            # 该点在某个三角形顶点上
            if relation == 1:
                concentration[i, j, :] = original_concentration[c_id]
            # 该点在某个三角形内部（包括边界）
            elif relation == 2:
                approximate_concentration = 0.0
                for k in range(0, 3):
                    nod_k = nods[c_id, k]
                    f_k = original_concentration[nod_k]
                    a_k = a_arr[c_id, k]
                    b_k = b_arr[c_id, k]
                    c_k = c_arr[c_id, k]
                    approximate_concentration = approximate_concentration + f_k * (a_k + b_k * x + c_k * y)
                concentration[i, j, :] = approximate_concentration
    return concentration


'''
        nano_original_concentration：纳米空间原始浓度
        lower_r/upper_r：根据该半径计算在纳米空间中的平均浓度
        xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
        radius_list：开放空间r方向上0~300的坐标
        relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
        a_arr,b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
        nods：保存纳米空间每个三角形三个点的编号
        '''


# 计算某个r<300时，该半径的平均值
def same_radius_avg(nano_original_concentration, r, xy_index, radius_list, relations, a_arr, b_arr, c_arr, nods):
    # 找到该半径在开放空间r的0~300范围内的下标
    indexes = np.where(radius_list == r)
    a = indexes[0][0]
    r_xy_index = xy_index[a]
    con_sum = 0.0
    count = 0
    for q in range(len(r_xy_index)):
        # 纳米空间中没有与当前xy对应的了
        if r_xy_index[q][0] == -1:
            break
        # r_xy_index[q][0]，r_xy_index[q][1]都是数组中的顺序，需要-300即为他们的xy值
        x = r_xy_index[q][0] - 200
        y = r_xy_index[q][1] - 200
        # 根据xy在纳米空间中插值计算该点的值
        concentration = cal_nano_concentration(x, y, nano_original_concentration, relations, a_arr, b_arr, c_arr, nods)
        # 累加
        con_sum = con_sum + concentration
        count = count + 1
    # 求平均值
    con_avg = con_sum / count
    return con_avg


# 使用哈希表（字典）来加速查找过程
def find_point_index_v2(coordinates_dict, target_coordinates):
    if target_coordinates in coordinates_dict:
        index = coordinates_dict[target_coordinates]
        return index
    else:
        print("Target coordinates not found in the array.")


'''
                   x_i, y_j：当前点三维坐标下xy
                   nano_original_concentration：纳米空间原始浓度
                   relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
                   a_arr, b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
                   nods：保存纳米空间每个三角形三个点的编号
'''


def cal_nano_concentration(x_i, y_j, nano_original_concentration, relations, a_arr, b_arr, c_arr, nods):
    # x_i, y_j为实际坐标，i、j为数组索引
    i = int(x_i) + 200
    j = int(y_j) + 200
    # 该点在三角形点上还是内部（包括边界）
    relation = relations[i, j, 0]
    # 三角形编号或点编号
    c_id = relations[i, j, 1]
    concentration = 0
    # 该点在某个三角形顶点上
    if relation == 1:
        # 直接取纳米空间网格点对应浓度
        concentration = nano_original_concentration[c_id]
    # 该点在某个三角形内部（包括边界）
    elif relation == 2:
        # 插值计算浓度
        approximate_concentration = 0.0
        for k in range(0, 3):
            nod_k = nods[c_id, k]
            f_k = nano_original_concentration[nod_k]
            a_k = a_arr[c_id, k]
            b_k = b_arr[c_id, k]
            c_k = c_arr[c_id, k]
            approximate_concentration = approximate_concentration + f_k * (a_k + b_k * x_i + c_k * y_j)
        concentration = approximate_concentration
    return concentration


'''
                    lower_r, upper_r, lower_z, upper_z：包含当前点网格的四个点的r与z
                    radius, height：当前点半径和高
                    open_original_concentration：开放空间原始浓度
                    coordinates_dict：哈希表（字典）用于存储坐标值与索引的映射
                    nano_original_concentration：纳米空间原始浓度
                    xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
                    radius_list：开放空间r方向上0~300的坐标
                    relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
                    a_arr, b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
                    nods：保存纳米空间每个三角形三个点的编号
'''


def interpolation_calculation(lower_r, upper_r, lower_z, upper_z, radius, height, open_original_concentration,
                              coordinates_dict, nano_original_concentration, xy_index, radius_list, relations, a_arr,
                              b_arr, c_arr, nods):
    # 使用双线性插值法
    r1, z1 = lower_r, lower_z  # 包含计算点的网格的左下点
    r2, z2 = upper_r, lower_z  # 包含计算点的网格的右下点
    r3, z3 = upper_r, upper_z  # 包含计算点的网格的右上点
    r4, z4 = lower_r, upper_z  # 包含计算点的网格的左上点
    alpha = (radius - lower_r) / (upper_r - lower_r)
    beta = (height - lower_z) / (upper_z - lower_z)
    # 该网格底部为开放空间与纳米空间上侧交界处（开放网格中没有这些浓度）
    if lower_z == 7.5 and radius < 200:
        # 找到右上点和左上点rz所对应的开放空间的点的下标（在字典里查找）
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        # v1, v2, v3, v4为浓度
        # same_radius_avg：因为该网格在交界处而开放空间中没有该点的浓度，所以在纳米空间中计算当前半径在纳米空间中的平均值，即为开放空间中该点的值
        '''
        nano_original_concentration：纳米空间原始浓度
        lower_r/upper_r：根据该半径计算在纳米空间中的平均浓度
        xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
        radius_list：开放空间r方向上0~300的坐标
        relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
        a_arr,b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
        nods：保存纳米空间每个三角形三个点的编号
        '''
        # 右上点和左上点浓度直接取开放空间网格点浓度，右下点和左下点浓度根据半径在纳米空间内求平均值
        v1, v2, v3, v4 = same_radius_avg(nano_original_concentration, lower_r, xy_index, radius_list, relations, a_arr,
                                         b_arr, c_arr, nods), same_radius_avg(
            nano_original_concentration, upper_r, xy_index, radius_list, relations, a_arr, b_arr, c_arr, nods), \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 该网格左侧在纳米空间内，右侧在交界处（200）
    elif radius == 200 and height <= 7.5:
        # index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        # 找到右上点rz所对应的开放空间的点的下标（在字典里查找）
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        # vl为该结构网格左侧两点的值（均相同）
        vl = same_radius_avg(nano_original_concentration, lower_r, xy_index, radius_list, relations, a_arr, b_arr,
                             c_arr, nods)
        # v1, v2, v3, v4 = vv, open_original_concentration[index2], open_original_concentration[index3], vv
        # 右下点浓度取纳米空间r=300处的平均值
        # 右上点浓度直接取开放空间网格点浓度，右下点、左下点和左上点浓度根据半径在纳米空间内求平均值
        v1, v2, v3, v4 = vl, same_radius_avg(nano_original_concentration, upper_r, xy_index, radius_list, relations,
                                             a_arr, b_arr,
                                             c_arr, nods), open_original_concentration[index3], vl

    # 该网格左下在上侧交界处
    elif radius == 200 and height <= 13.0:
        # 找到右下、右上、左上点rz所对应的开放空间的点的下标（在字典里查找）
        index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        # 右下、右上、左上点浓度直接取开放空间网格点浓度，左下点浓度根据半径在纳米空间内求平均值
        v1, v2, v3, v4 = same_radius_avg(nano_original_concentration, lower_r, xy_index, radius_list, relations, a_arr,
                                         b_arr, c_arr, nods), open_original_concentration[
                             index2], \
                         open_original_concentration[index3], open_original_concentration[index4]
    # 四个点的浓度都在开放空间文件中
    else:
        # 找到四个点rz所对应的开放空间的点的下标（在字典里查找）
        index1 = find_point_index_v2(coordinates_dict, (z1, r1))
        index2 = find_point_index_v2(coordinates_dict, (z2, r2))
        index3 = find_point_index_v2(coordinates_dict, (z3, r3))
        index4 = find_point_index_v2(coordinates_dict, (z4, r4))
        # 四个点直接取开放空间网格点浓度
        v1, v2, v3, v4 = open_original_concentration[index1], open_original_concentration[index2], \
                         open_original_concentration[index3], \
                         open_original_concentration[index4]
    # 插值计算浓度
    interpolated_value = (1 - alpha) * (1 - beta) * v1 + alpha * (1 - beta) * v2 + alpha * beta * v3 + (
            1 - alpha) * beta * v4
    return interpolated_value


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


def process_concentration(nano_original_concentration, open_original_concentration, c_val, position_list_i, grids_zr,
                          xy_index,
                          radius_list, coordinates_dict, relations, a_arr, b_arr, c_arr, nods):
    print("process_concentration 开始")
    dis_r = position_list_i[0]  # 0/300/500，当前观察点
    dis_z = position_list_i[1]  # 0/300/500，当前观察点
    # 以当前观察点为基准，向四周根据不同距离，布浓度矩阵点。half_length=500
    x_arr, y_arr, z_arr, process_points = generate_interval(dis_r, dis_z, half_length)
    d1, d2, d3 = len(x_arr), len(y_arr), len(z_arr)  # 都是83
    # x和y平方的数组
    x_square, y_square = np.square(x_arr), np.square(y_arr)
    # xyz绝对值的数组，因为对称性，所以只需正值即可
    x_abs, y_abs, z_abs = np.abs(x_arr), np.abs(y_arr), np.abs(z_arr)
    # 存储最终要生成的浓度矩阵
    concentration = np.full((d1, d2, d3), fill_value=c_val)
    # x
    for i in range(d1):
        # 用于生成半径
        x2 = x_square[i]
        # grids_zr中的下标
        a = x_abs[i]
        # 用于计算在纳米空间的点
        x_i = x_arr[i]
        # y
        for j in range(d2):
            # 用于生成半径
            y2 = y_square[j]
            # grids_zr中的下标
            b = y_abs[j]
            # 用于计算在纳米空间的点
            y_j = y_arr[j]
            # z
            for k in range(d3):
                # 高度z
                height = z_abs[k]
                # 半径r
                radius = np.sqrt(x2 + y2)
                # grids_zr中的下标
                c = z_abs[k]
                # 若该点在纳米空间外
                if height >= 7.5 or radius >= 200:
                    # 取出将其包含在内的网格的rz
                    lower_r, upper_r, lower_z, upper_z = grids_zr[np.int32(a), np.int32(b), np.int32(c)]
                    # 插值计算开放空间待插入点浓度
                    '''
                    lower_r, upper_r, lower_z, upper_z：包含当前点网格的四个点的r与z
                    radius, height：当前点半径和高
                    open_original_concentration：开放空间原始浓度
                    coordinates_dict：哈希表（字典）用于存储坐标值与索引的映射
                    nano_original_concentration：纳米空间原始浓度
                    xy_index：该文件保存浓度矩阵中，布在纳米空间内的点，与在纳米空间中对应的网格点的i和j（顺序）
                    radius_list：开放空间r方向上0~300的坐标
                    relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
                    a_arr, b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
                    nods：保存纳米空间每个三角形三个点的编号
                    '''
                    concentration[i, j, k] = interpolation_calculation(lower_r, upper_r, lower_z, upper_z,
                                                                       radius,
                                                                       height,
                                                                       open_original_concentration,
                                                                       coordinates_dict, nano_original_concentration,
                                                                       xy_index,
                                                                       radius_list, relations, a_arr, b_arr, c_arr,
                                                                       nods)
                # 插值计算纳米空间待插入点浓度
                else:
                    '''
                    x_i, y_j：当前点三维坐标下xy
                    nano_original_concentration：纳米空间原始浓度
                    relations：保存浓度矩阵中，布在纳米空间内的点，所在的非结构网格编号或网格点编号
                    a_arr, b_arr, c_arr：对应公式中的ai bi ci，用于非结构网格中的插值计算
                    nods：保存纳米空间每个三角形三个点的编号
                    '''
                    concentration[i, j, k] = cal_nano_concentration(x_i, y_j, nano_original_concentration, relations,
                                                                    a_arr, b_arr, c_arr, nods)
    # 返回浓度矩阵与中心点位置
    return concentration, len(process_points)
