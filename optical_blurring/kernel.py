import numpy as np
from scipy.ndimage import convolve1d
from ob_parameters import *


# 高斯函数
def gaussian_function(w, sigma):
    return np.power((2 * np.pi * sigma), -1 / 2) * np.exp(-np.square(w) / (2 * sigma))


# 生成卷积核
def generate_kernel():
    # sigma_xy = np.square(400) / (8 * np.log(2))
    sigma = sigma_xy
    # x、y方向上，xy_count=41
    xy_parameter = np.array([gaussian_function(i, sigma) for i in range(xy_count + 1)])
    # sigma_z = np.square(800) / (8 * np.log(2))
    sigma = sigma_z
    # z方向上，z_count=41
    z_parameter = np.array([gaussian_function(i, sigma) for i in range(z_count + 1)])
    # xy方向拼接成83
    xy_parameter = np.concatenate((xy_parameter[xy_count + 1:0:-1], xy_parameter))
    xy_parameter /= np.sum(xy_parameter)  # 归一化
    # z方向拼接成83
    z_parameter = np.concatenate((z_parameter[z_count + 1:0:-1], z_parameter))
    z_parameter /= np.sum(z_parameter)  # 归一化
    # 得到卷积核矩阵大小为83*83*83，因为浓度矩阵大小为83*83*83
    kernel = np.empty((len(xy_parameter), len(xy_parameter), len(z_parameter)))
    for i in range(len(xy_parameter)):
        for j in range(len(xy_parameter)):
            for k in range(len(z_parameter)):
                kernel[i, j, k] = xy_parameter[i] * xy_parameter[j] * z_parameter[k]
    np.set_printoptions(precision=20)
    return kernel


'''
original_matrix：浓度矩阵 83*83*83
kernel：卷积核矩阵 165*165*165，即generate_kernel()生成的矩阵
xy_c_val：0.0 用于纳米钙火花 /4.081632653061224858e-03 用于钙火花
z_c_val：0.0 用于纳米钙火花 /4.081632653061224858e-03 用于钙火花，
'''


# 三维卷积
def convolve3d(original_matrix, kernel, xy_c_val, z_c_val):
    # xy上的填充值
    c_val = xy_c_val
    # 维度，从0开始
    i = 0
    print("convolve3d开始")
    # 从xyz三个维度进行卷积
    for kernel_i in kernel:
        # 在z上卷积时，更换填充值
        if i == 2:
            c_val = z_c_val
        # 在第一个维度上进行一维卷积操作
        '''
        from scipy.ndimage import convolve1d
        original_matrix：浓度矩阵 83*83*83
        kernel_i：卷积核的一维数组，与 original_matrix 在指定的轴 (axis=i) 上进行卷积操作。即第i个维度的卷积核矩阵，大小为165 
        axis=i：指定在哪个轴上对浓度矩阵进行一维卷积
        mode='constant'：在浓度矩阵的边界外部使用常数值进行填充
        cval=c_val：填充浓度矩阵的边界
        '''
        original_matrix = convolve1d(original_matrix, kernel_i, axis=i, mode='constant', cval=c_val)
        # 维度加一
        i = i + 1
    return original_matrix


def main():
    kernel = generate_kernel()
    np.save("kernel/kernel_3D", kernel)
    # a=np.load("kernel/kernel_3D.npy")


if __name__ == "__main__":
    main()
