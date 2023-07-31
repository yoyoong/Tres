import numpy as np
import pandas as pd

import pywt

def getData():
    train_x = np.load('./A/train_x.npz')['data']  # (21547, 4, 15000)
    train_y = np.load('./A/train_y.npz')['data']  # （21547,)
    test_X = np.load('./A/test_x.npz')['data']  # (1000, 4, 15000)
    return train_x, train_y, test_X

if __name__ == '__main__':
    train_x, train_y, submit_X = getData()
    test_x = train_x[20000:train_x.shape[0]]
    train_x = train_x[0:20000]
    test_y = train_y[20000:train_y.shape[0]]
    train_y = train_y[0:20000]

    # train_x_y0 = train_x[train_y == 0]
    # train_x_y1 = train_x[train_y == 1]
    # train_x_y2 = train_x[train_y == 2]
    # train_x_y3 = train_x[train_y == 3]

    for n in range(test_x.shape[0]):
        rmean_array = np.zeros([4, 4], dtype=float)
        rtop1000mean_array = np.zeros([4, 4], dtype=float)
        rmax_array = np.zeros([4, 4], dtype=float)
        for i in range(train_x.shape[1]):
            train_x_ci = train_x[:, i, :]
            for j in range(4):
                train_x_ci_yj = train_x_ci[train_y == j]
                train_x_ci_yj_add = np.row_stack((train_x_ci_yj, train_x[n, i, :]))
                r_matrix = np.corrcoef(train_x_ci_yj_add)[:, -1]
                rmean_array[i][j] = np.mean(r_matrix)
                rmax_array[i][j] = np.sort(r_matrix)[-2]
                rtop1000mean_array[i][j] = np.mean(abs(np.sort(-r_matrix))[0:1000])
        rmaxmean = np.mean(rmax_array, axis=0)
        rtop1000meanmean = np.mean(rtop1000mean_array, axis=0)
        rmaxmean_max = np.argmax(rmaxmean)
        rtop1000meanmean_max = np.argmax(rtop1000meanmean)
        print("rmaxmean_max:", rmaxmean_max, "rtop1000meanmean_max:", rtop1000meanmean_max, "y:", test_y[n])



    # for i in range(train_x.shape[1]):
    #     train_x_df = pd.DataFrame(train_x[:, i, :].reshape(train_x.shape[0], train_x.shape[2]))
    #     print(train_x_df.T.describe(include='all'))





    wavelet = 'db2'
    for i in range(train_x.shape[0]):
        for j in range(train_x.shape[1]):
            signal = train_x[i][j]
            coeffs = pywt.wavedec(train_x[i][j], wavelet)
            print(coeffs)

    print(coeffs)

