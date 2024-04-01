from nano_parameters import *
import numpy as np
from tool.cal_bcnl import *


def set_cag_0(cag_con):
    cag = np.copy(cag_con)
    points_index = np.loadtxt("r_gt_100.txt", dtype=int)
    for i in points_index:
        cag[i] = 0.0
    return cag


def cal_dye_and_buffers(f, caf, cag, cab1, cab2, cab3, cab4):
    j_fdye = np.zeros(NP)
    j_gdye = np.zeros(NP)
    j_1 = np.zeros(NP)
    j_2 = np.zeros(NP)
    j_3 = np.zeros(NP)
    j_4 = np.zeros(NP)
    new_cag = np.zeros(NP)
    new_cab1 = np.zeros(NP)
    new_cab2 = np.zeros(NP)
    new_cab3 = np.zeros(NP)
    new_cab4 = np.zeros(NP)
    for i in range(0, NP):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
        j_gdye[i] = -K_GCaMP6f_PLUS * f[i] * (GCaMP6f_T - cag[i]) + K_GCaMP6f_MINUS * cag[i]
        j_1[i] = -K_Calmodulin_PLUS * f[i] * (Calmodulin_T - cab1[i]) + K_Calmodulin_MINUS * cab1[i]
        j_2[i] = -K_TroponinC_PLUS * f[i] * (TroponinC_T - cab2[i]) + K_TroponinC_MINUS * cab2[i]
        j_3[i] = -K_SR_PLUS * f[i] * (SR_T - cab3[i]) + K_SR_MINUS * cab3[i]
        j_4[i] = -K_SL_PLUS_NANO * f[i] * (SL_T_NANO - cab4[i]) + K_SL_MINUS_NANO * cab4[i]

    for i in range(0, NP):
        new_cag[i] = (-j_gdye[i]) * DT + cag[i]
        new_cab1[i] = (-j_1[i]) * DT + cab1[i]
        new_cab2[i] = (-j_2[i]) * DT + cab2[i]
        new_cab3[i] = (-j_3[i]) * DT + cab3[i]
        new_cab4[i] = (-j_4[i]) * DT + cab4[i]
    return j_fdye, j_gdye, j_1, j_2, j_3, j_4, new_cag, new_cab1, new_cab2, new_cab3, new_cab4


'''
        k_ryr：ryr通道处钙离子释放的扩散系数
        nano_f, nano_caf, nano_cag, nano_cab3, nano_cab4：当前步数纳米空间ca、caf、cag、SR membrane、SL membrane浓度
        ca_jsr：jSR中的钙离子浓度
        c_ca_out：边界浓度
        bcnl_elements：和三角形有关的参数
'''


# 计算Ca
def nano_calculation_f(k_ryr, f, caf, cag, cab1, cab2, cab3, cab4, ca_jsr, c_ca_out, bcnl_elements):
    # 计算染料和缓冲物
    j_fdye, j_gdye, j_1, j_2, j_3, j_4, new_cag, new_cab1, new_cab2, new_cab3, new_cab4 = cal_dye_and_buffers(f, caf,
                                                                                                              cag, cab1,
                                                                                                              cab2,
                                                                                                              cab3,
                                                                                                              cab4)
    j_gdye = set_cag_0(j_gdye)
    new_cag = set_cag_0(new_cag)
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
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = bcnl_elements
    # 保存此次迭代结果
    last = np.zeros(NP)
    # 临时保存结果
    temp = np.zeros(NP)
    # 迭代十次
    for j in range(0, 10):
        for i in range(0, NP):
            item_1 = (j_fdye[i] + j_gdye[i] + j_1[i] + j_2[i] + j_3[i] + j_4[i]) * control_area[i]
            # S*{J_Gdye+J_Fdye+J_buffers }
            item_2, item_3, item_4, item_5, item_6, item_7 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            # nmax=9
            for k in range(0, nmax):
                # 得到控制面积内的三角形编号
                near_triangle_index = near_triangle[i, k]
                # 若有三角形
                if near_triangle_index != -1:
                    # 该点i在该三角形（near_triangle_index）中的编号
                    point_index = index_in_triangle[i, k]
                    # 另外两个点的在三角形中编号
                    point2 = (point_index + 1) % 3
                    point3 = (point_index + 2) % 3
                    # 该点i所对应的边的两点的编号
                    nod2 = nod[near_triangle_index, point2]
                    nod3 = nod[near_triangle_index, point3]
                    '''
                    in_boundary：若有入流边，则存储边界对应的的两个点
                    out_boundary：若有出流边，则存储边界对应的的两个点
                    inner：否则三个点都存在里面
                    '''
                    in_boundary, out_boundary, inner = judge_point(near_triangle_index, nod, npoch)

                    # 求导部分的项没有边界边（当前点对应的边不为边界）
                    if npoch[nod2] == B_INNER or npoch[nod3] == B_INNER:
                        # 若第一次迭代：求导、出流边界、入流边界中使用的除所求点之外的其余两个顶点的值使用的n∆t时刻的值，即f_n
                        f2i, f3i = f[nod2], f[nod3]
                        # 若不是第一次迭代，求导中使用的除所求点之外的其余两个顶点的值使用上一时刻的值与上一次迭代值得平均值
                        if j != 0:
                            # 学长只取上一时刻的值
                            f2i, f3i = last[nod2], last[nod3]
                            # f2i, f3i = (f[nod2] + last[nod2]) / 2.0, (f[nod3] + last[nod3]) / 2.0
                        # D_ca*∑{[(f_2i*b_2i+f_3i*b_3i)*(e_x) ⃗(n_xi) ⃗+(f_2i*c_2i+f_3i*c_3i)*(e_y) ⃗(n_yi) ⃗ ]∙L_i}
                        # nix*l niy*l之前已经计算得出
                        item_2 = item_2 + D_CA * (
                                (f2i * b_arr[near_triangle_index, point2] + f3i * b_arr[near_triangle_index, point3]) *
                                nix_multiply_l[near_triangle_index, point_index] + (
                                        f2i * c_arr[near_triangle_index, point2] + f3i * c_arr[
                                    near_triangle_index, point3]) * niy_multiply_l[
                                    near_triangle_index, point_index])
                        # D_ca*∑{[b_1i*(e_x) ⃗(n_xi) ⃗+c_1i*(e_y) ⃗(n_yi) ⃗ ]∙L_i}
                        item_3 = item_3 + D_CA * (b_arr[near_triangle_index, point_index] * nix_multiply_l[
                            near_triangle_index, point_index] + c_arr[near_triangle_index, point_index] *
                                                  niy_multiply_l[
                                                      near_triangle_index, point_index])
                    # 若不是第一次迭代，出流边界、入流边界中使用的除所求点之外的其余两个顶点的值使用上一次迭代得到的值，即f_(p-1)
                    f2jk, f3jk = f[nod2], f[nod3]
                    if j != 0:
                        f2jk, f3jk = last[nod2], last[nod3]
                    fjk = (f2jk + f3jk) / 3.0
                    # 若该三角形为入流边界，并且通道未关闭
                    if len(in_boundary) == 2 and k_ryr != 0:
                        # 计算边界边长
                        length_j = cal_length(in_boundary, grid)
                        # K_ryr*∑([〖Ca〗^(2+)]_jSR-(f_2j+f_3j)/3)*L_j
                        item_4 = item_4 + k_ryr * (ca_jsr - fjk) * length_j
                        # -(K_ryr*∑L_j)/3
                        item_5 = item_5 - k_ryr * length_j / 3.0
                    # 若该三角形为出流边界
                    elif len(out_boundary) == 2:
                        # 计算边界边长
                        length_k = cal_length(out_boundary, grid)
                        # (D_ca/∆r)*∑([〖Ca〗^(2+)]_out-(f_2k+f_3k)/3)*L_k
                        # item_6 = item_6 + (D_CA / Delta_r) * (c_ca_out - fjk) * length_k
                        item_6 = item_6 + (D_CA_OUT / Delta_r) * (c_ca_out - fjk) * length_k
                        # -(D_ca/∆r)*∑L_k)/3
                        # item_7 = item_7 - (D_CA / Delta_r) * length_k / 3.0
                        item_7 = item_7 - (D_CA_OUT / Delta_r) * length_k / 3.0
            # 当前点(n+1)Δt时刻的值，存入temp中，不会影响本次迭代其他点的计算
            temp[i] = ((item_2 + item_1 + item_4 + item_6) * DT + control_area[i] * f[i]) / (
                    control_area[i] - (item_3 + item_5 + item_7) * DT)
        # 此次迭代的值保存到last数组
        last = np.copy(temp)
    # 这一步的值保存到new_f
    new_f = np.copy(last)
    return new_f, new_cag, new_cab1, new_cab2, new_cab3, new_cab4


'''
        nano_f, nano_caf：当前步数纳米空间ca、caf浓度
        c_caf_out：边界浓度
        bcnl_elements：和三角形有关的参数
'''


# CaF
def nano_calculation_g2(f, caf, c_caf_out, bcnl_elements):
    # 染料流量
    j_fdye = np.zeros(NP)
    for i in range(0, NP):
        j_fdye[i] = -K_F3_PLUS * f[i] * (F3_T - caf[i]) + K_F3_MINUS * caf[i]
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
    single_area, control_area, near_triangle, index_in_triangle, nix_multiply_l, niy_multiply_l, a_arr, b_arr, c_arr, nmax, total_area = bcnl_elements
    # 保存此次迭代结果
    last = np.zeros(NP)
    # 临时保存结果
    temp = np.zeros(NP)
    # 迭代十次
    for j in range(0, 10):
        # 遍历每个点
        for i in range(0, NP):
            # S*(-J_Fdye)
            item_1 = -j_fdye[i] * control_area[i]
            item_2, item_3, item_4, item_5 = 0.0, 0.0, 0.0, 0.0
            for k in range(0, nmax):
                # 得到控制面积内的三角形编号
                near_triangle_index = near_triangle[i, k]
                # 若有三角形
                if near_triangle_index != -1:
                    # 该点在该三角形（near_triangle_index）中的编号
                    point_index = index_in_triangle[i, k]
                    # 另外两个点的在三角形中编号
                    point2 = (point_index + 1) % 3
                    point3 = (point_index + 2) % 3
                    # 该点i所对应的边的两点的编号
                    nod2 = nod[near_triangle_index, point2]
                    nod3 = nod[near_triangle_index, point3]
                    '''
                        in_boundary：若有入流边，则存储边界对应的的两个点
                        out_boundary：若有出流边，则存储边界对应的的两个点
                        inner：否则三个点都存在里面
                    '''
                    in_boundary, out_boundary, inner = judge_point(near_triangle_index, nod, npoch)
                    # 求导部分的项没有边界边（当前点对应的边不为边界）
                    if npoch[nod2] == B_INNER or npoch[nod3] == B_INNER:
                        # 若第一次迭代：求导、出流边界、入流边界中使用的除所求点之外的其余两个顶点的值使用的n∆t时刻的值，即g_n
                        g2i, g3i = caf[nod2], caf[nod3]
                        # 若不是第一次迭代，求导中使用的除所求点之外的其余两个顶点的值使用上一时刻的值与上一次迭代值得平均值
                        if j != 0:
                            # 学长只取上一时刻的值
                            g2i, g3i = last[nod2], last[nod3]
                            # g2i, g3i = (caf[nod2] + last[nod2]) / 2.0, (caf[nod3] + last[nod3]) / 2.0
                        # D_caf*∑{[(g_2i*b_2i+g_3i*b_3i)*(e_x) ⃗(n_xi) ⃗+(g_2i*c_2i+g_3i*c_3i)*(e_y) ⃗(n_yi) ⃗ ]∙L_i}
                        # nix*l niy*l之前已经计算得出
                        item_2 = item_2 + D_CAF * ((g2i * b_arr[near_triangle_index, point2] + g3i * b_arr[
                            near_triangle_index, point3]) * nix_multiply_l[near_triangle_index, point_index] + (
                                                           g2i * c_arr[near_triangle_index, point2] + g3i * c_arr[
                                                       near_triangle_index, point3]) * niy_multiply_l[
                                                       near_triangle_index, point_index])
                        # D_caf*∑{[b_1i*(e_x) ⃗(n_xi) ⃗+c_1i*(e_y) ⃗(n_yi) ⃗ ]∙L_i}
                        item_3 = item_3 + D_CAF * (b_arr[near_triangle_index, point_index] * nix_multiply_l[
                            near_triangle_index, point_index] + c_arr[near_triangle_index, point_index] *
                                                   niy_multiply_l[near_triangle_index, point_index])
                    # 若不是第一次迭代，出流边界、入流边界中使用的除所求点之外的其余两个顶点的值使用上一次迭代得到的值，即g_(p-1)
                    g2k, g3k = caf[nod2], caf[nod3]
                    if j != 0:
                        g2k, g3k = last[nod2], last[nod3]
                    gk = (g2k + g3k) / 3.0
                    # 若该三角形为出流边界
                    if len(out_boundary) == 2:
                        # 计算边界边长
                        length_k = cal_length(out_boundary, grid)
                        # (D_caf/∆r)*∑([〖Ca〗^(2+)]_out-(g_2k+g_3k)/3)*L_k
                        # item_4 = item_4 + (D_CAF / Delta_r) * (c_caf_out - gk) * length_k
                        item_4 = item_4 + (D_CAF_OUT / Delta_r) * (c_caf_out - gk) * length_k
                        # -(D_caf/∆r)*∑L_k)/3
                        # item_5 = item_5 - (D_CAF / Delta_r) * length_k / 3.0
                        item_5 = item_5 - (D_CAF_OUT / Delta_r) * length_k / 3.0
            # 当前点(n+1)Δt时刻的值，存入temp中，不会影响本次迭代其他点的计算
            temp[i] = ((item_2 + item_1 + item_4) * DT + control_area[i] * caf[i]) / (
                    control_area[i] - (item_3 + item_5) * DT)
        # 此次迭代的值保存到last数组
        last = np.copy(temp)
    # 这一步的值保存到new_caf
    new_caf = np.copy(last)
    return new_caf
