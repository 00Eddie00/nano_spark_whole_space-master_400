import numpy as np
from shapely.geometry import Point, Polygon

from tool.point_scatter_util import point_scatter


def nano_judge_relation():
    relations = np.full((401, 401, 2), dtype=int, fill_value=-1)
    # 创建点坐标数组
    nods = np.loadtxt("../config/nano/4RYRnod.dat", dtype=int) - 1
    points = np.loadtxt("../config/nano/4RYRgridt.dat", dtype=np.float64)
    triangle_items = []
    for i in range(len(nods)):
        triangle1_vertices = []
        for j in range(3):
            nod_j = nods[i, j]
            x, y = points[nod_j, 0], points[nod_j, 1]
            triangle1_vertices.append((x, y))
        triangle_i = Polygon(triangle1_vertices)
        triangle_items.append(triangle_i)
    # 创建三角形集合对象
    print("三角形集合对象创建完成")
    for i in range(401):
        print(f"i={i}")
        for j in range(401):
            x = i - 200
            y = j - 200
            if np.square(x) + np.square(y) <= 40000:
                # 创建要插入的新点对象
                new_point = Point(x, y)
                # 检查新点是否在三角形集合中的任何一个三角形内部
                for p, triangle in enumerate(triangle_items):
                    # 检查新点是否在三角形内部，包括边界（顶点）
                    if triangle.contains(new_point) or triangle.touches(new_point):
                        triangle_index = p  # 记录包含新点的三角形索引
                        relations[i, j, 0] = 2
                        relations[i, j, 1] = triangle_index
                        # 检查该点是否在顶点上
                        for k in range(3):
                            nod_k = nods[triangle_index, k]
                            xx, yy = points[nod_k, 0], points[nod_k, 1]
                            if x == xx and y == yy:
                                relations[i, j, 0] = 1
                                relations[i, j, 1] = nod_k
                                break
                        break
    return relations


# 生成开放空间的点的相邻点
def nano_relations_refine():
    relations = np.load("relation/relations.npy")
    refined_relations = np.copy(relations)
    print("修改一些点")
    for i in range(401):
        for j in range(401):
            x = i - 200
            y = j - 200
            r_square = np.square(x) + np.square(y)
            if r_square <= 40000 and relations[i, j, 0] == -1:
                r = np.sqrt(r_square)
                if r > 199:
                    refined_relations[i, j, 0] = 1
                    refined_relations[i, j, 1] = 84
                elif r < 199:
                    refined_relations[i, j, 0] = 1
                    if x == 15 and y == 15:
                        refined_relations[i, j, 1] = 0
                    elif x == 15 and y == -15:
                        refined_relations[i, j, 1] = 21
                    elif x == -15 and y == -15:
                        refined_relations[i, j, 1] = 42
                    elif x == -15 and y == 15:
                        refined_relations[i, j, 1] = 63
    np.save("refined_relations", refined_relations)



def cal_nano_points():
    new_open_r = np.flipud(point_scatter(195, 0, 5, k=2, positive=False))
    radius_list = np.concatenate((new_open_r, (
        200,)))
    xy_index = np.full((len(radius_list), 50, 2), fill_value=-1)
    xy_len_index = np.zeros(len(radius_list), dtype=int)
    for i in range(401):
        x = i - 200
        x2 = np.square(x)
        for j in range(401):
            y = j - 200
            y2 = np.square(y)
            r = np.sqrt(x2 + y2)
            # 找到需要布的点和纳米空间中的点对应点的下标
            indexes = np.where(radius_list == r)
            if len(indexes[0]) > 0:
                # 对应点的下标
                a = indexes[0][0]
                # 记录该下标存了几个
                b = xy_len_index[a]
                # 存储对应的在x、y中的顺序下标
                xy_index[a][b][0], xy_index[a][b][1] = i, j
                # 存储个数加一
                xy_len_index[a] = b + 1
    np.save(f"../optical_blurring/xy_list/xy_index", xy_index)


def main():
    # relations = nano_judge_relation()
    # np.save("relation/relations", relations)
    # nano_relations_refine()
    cal_nano_points()


if __name__ == "__main__":
    main()
