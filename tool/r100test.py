import numpy as np


def gt100():
    cor = np.loadtxt("../config/nano/4RYRgridt.dat")
    l = 0
    temp = np.empty(1131)
    for j in range(len(cor)):
        c = cor[j]
        r = np.sqrt(c[0] * c[0] + c[1] * c[1])
        if r >= 100.0:
            temp[l] = j
            print(j, r, c[0], c[1])
            l = l + 1
    # np.savetxt("r_gt_100.txt", temp, fmt='%d')
    print(l)


def main():
    gt100()


if __name__ == "__main__":
    main()
