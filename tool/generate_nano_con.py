import numpy as np


# dis_r, dis_z = 300, 300
# half_length=500
def generate_interval(dis_r, dis_z, half_length):
    # 设置间隔，每隔一定间隔取一个点
    interval_arr = [5, 10, 20, 25, 50]
    end_value = 0.0
    # 每隔100nm，间隔不同
    interval = half_length // len(interval_arr)
    process_points = np.empty(0)
    for i in range(len(interval_arr)):
        start_value = end_value
        end_value = end_value + interval
        num_points = interval // interval_arr[i] + 1
        data_array = np.linspace(start_value, end_value, int(num_points))
        if i != len(interval_arr) - 1:
            data_array = data_array[:-1]
        process_points = np.concatenate((process_points, data_array))
    # 生成
    x_half_left = dis_r - np.flip(process_points)
    x_half_right = process_points + dis_r
    x_arr = np.concatenate((x_half_left[:-1], x_half_right))

    y_half_left = - np.flip(process_points)
    y_half_right = process_points
    y_arr = np.concatenate((y_half_left[:-1], y_half_right))

    z_half_left = dis_z - np.flip(process_points)
    z_half_right = process_points + dis_z
    z_arr = np.concatenate((z_half_left[:-1], z_half_right))
    return x_arr, y_arr, z_arr, process_points

