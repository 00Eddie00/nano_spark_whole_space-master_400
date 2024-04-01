import numpy as np
import os


def cal_ctrl_v():
    grid_coordinates = np.loadtxt("../../config/open/open_grid_coordinates.csv", delimiter=",")
    neighbors = np.loadtxt("../../config/open/open_neighbor.csv", int, delimiter=",")
    NP = len(grid_coordinates)
    ctrl_v = np.zeros(NP)
    for i in range(NP):
        z2 = grid_coordinates[i][0]
        r2 = grid_coordinates[i][1]
        neighbor_up, neighbor_down, neighbor_right, neighbor_left = neighbors[i][0], neighbors[i][1], neighbors[i][2], \
                                                                    neighbors[i][3]  # 上下外内
        # 没有上顶面
        if neighbor_up == -1:
            # 则一定有下顶面
            z1 = grid_coordinates[neighbor_down, 0]
            h = z2 - z1
        # 没有下顶面
        elif neighbor_down == -1:
            # 则一定有上顶面
            z3 = grid_coordinates[neighbor_up, 0]
            h = z3 - z2
        else:
            z3 = grid_coordinates[neighbor_up, 0]
            z1 = grid_coordinates[neighbor_down, 0]
            h = z3 - z1

        # 没有外侧面
        if neighbor_right == -1:
            # 则一定有内侧面
            r1 = grid_coordinates[neighbor_left, 1]
            r3 = r2
        # 没有内侧面
        elif neighbor_left == -1:
            # 则一定有外侧面
            r3 = grid_coordinates[neighbor_right, 1]
            r1 = r2
        else:
            r1 = grid_coordinates[neighbor_left, 1]
            r3 = grid_coordinates[neighbor_right, 1]
        r1_5 = (r2 + r1) / 2
        r2_5 = (r2 + r3) / 2
        # h为整个高度，计算总体积时需要除以二，ctrl_v
        # ctrl_v[i] = np.pi * (r2_5 ** 2 - r1_5 ** 2) * h
        # h为整个高度除以二（即控制体积），计算总体积时不需要除以二，ctrl_v2
        # 其实是一样的。。。
        h = h / 2
        ctrl_v[i] = np.pi * (np.square(r2_5) - np.square(r1_5)) * h
    return ctrl_v



def main():
    np.save('ctrl_v2', cal_ctrl_v())

if __name__ == "__main__":
    main()
