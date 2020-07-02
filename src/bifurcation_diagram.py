import numpy as np
import pickle
import h5py
import os
import io


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
    if os.path.isfile('graph_list'):
        file = io.open('graph_list', 'rb')
        graph = pickle.load(file)
        file.close()
    else:
        graph = []
    A_step = 0.0001
    for i in range(100):
        velocities = []
        A = (i+1)*A_step + 0.04
        group = "group{:.4f}".format(A)
        with h5py.File('../data/run5/data.h5', 'r') as f:
            dset = f[group]
            print(A)
            for j in range(500):
                data = dset['dset' + str(j)]
                NY = int(len(data)/3)
                data = np.reshape(data, (3, NY))
                length = len(data[0])
                index = np.argmax(data[2][-length//4:-length//4 + 2001])
                start = length-length//4 + index
                y_vals = get_unique(data[2], start, 100)
                velocities = add_unique(velocities, y_vals)
        print(velocities)
        graph.append([round(A, 4), velocities])
    print(graph)
    file = io.open('graph_list', 'wb')
    pickle.dump(graph, file)
    file.close()


if __name__ == '__main__':
    main()
