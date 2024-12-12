
import numpy as np


base_dir = "/vast/home/yoonsoo/flastro/build/"

####################################################################
time_index = 10 * np.arange(31)

num_ranks = 36

####################################################################

data_stitched = None

for t_index in time_index:

    print(t_index)

    for i in np.arange(num_ranks):
        # print("i = ", i)

        file_name = base_dir + "output-%05d-2D-%d.raw" % (t_index, i)
        data_i = np.loadtxt(file_name).T

        if i == 0:
            data_stitched = data_i
        else:
            data_stitched = np.hstack((data_stitched, data_i))

    time, ix, iy, x, y, density, pressure, vx, vy, totalE, radE = data_stitched

    ind = np.lexsort((y, x))
    data_stitched = data_stitched[:, ind]

    time, ix, iy, x, y, density, pressure, vx, vy, totalE, radE = data_stitched

    np.savetxt("data_2D_%05d" % t_index, data_stitched)
