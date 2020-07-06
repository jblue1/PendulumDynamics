import numpy as np
import h5py
import matplotlib.pyplot as plt


def get_unique(arr, start, step):
    unique = []
    i = start
    while i < len(arr):
        if len(unique) < 1:
            unique.append(arr[i])
        else:
            diff = 0
            for j in unique:
                if np.abs(j - arr[i]) > 1e-6:
                    diff += 1
            if diff == len(unique):
                unique.append(arr[i])
        i += step
    return unique


def add_unique(list1, list2):
    unique = list1
    for i in list2:
        if len(unique) < 1:
            unique.append(i)
        else:
            diff = 0
            for j in unique:
                if np.abs(j - i) > 1e-6:
                    diff += 1
            if diff == len(unique):
                unique.append(i)
    return unique


def main():
    bif = []
    with h5py.File('../data/big/data.h5', 'r') as f:
        A = 0.0001
        for i in range(699):
            unique_thetas = []
            group = "group{:.4f}".format(A)
            dset = f[group]
            print(A)
            for j in range(500):
                data = dset['dset' + str(j)]
                NY = int(len(data)/3)
                data = np.reshape(data, (3, NY))
                y_vals = get_unique(data[2], 500, 1)
                unique_thetas = add_unique(unique_thetas, y_vals)
            bif.append([A, unique_thetas])
            A += 0.0001

    As = []
    vals = []
    for i in range(len(bif)):
        for j in range(len(bif[i][1])):
            As.append(bif[i][0])
            vals.append(bif[i][1][j])

    fig = plt.figure(figsize=(8, 5), dpi=300)
    ax = fig.add_subplot(111)
    ax.scatter(As, vals, s=0.5)
    ax.set_title('Bifurcation Diagram')
    ax.set_xlabel('A')
    ax.set_ylabel('Theta_dot')
    ax.set_xlim(0.015, 0.05)
    ax.set_ylim(-2, 2)

    plt.savefig('../data/imgs/Bifurcation_Diagram_Zoom.png')

    return


if __name__ == '__main__':
    main()
