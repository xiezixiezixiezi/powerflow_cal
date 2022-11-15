import numpy as np
import time
# from matplotlib import pyplot as plt
#
# plt.style.use('ggplot')

start_time = time.time()
# 形成导纳矩阵
g1 = -10
g2 = 0.01
Y = np.array([[2*(g1+g2), -g1, -g1],
              [-g1, 2*(g1+g2), -g1],
              [-g1, -g1, 2*(g1+g2)]])

print(Y)


# 置初值
U = np.array([1.0, 1.05, 1.0])
theta = np.array([0.0, 0.0, 0.0])

P = np.array([-2.8653, -0.6661])
Q = -1.2244

iteration = 0
T = 1

while T == 1:

    P_i = np.ones([1, 3])
    Q_i = np.ones([1, 3])

    # 计算潮流方程
    for i in range(3):
        P_i[0, i] = U[i] * (U[0] * Y[i, 0] * np.sin(theta[i] - theta[0])
                            + U[1] * Y[i, 1] * np.sin(theta[i] - theta[1])
                            + U[2] * Y[i, 2] * np.sin(theta[i] - theta[2]))
        Q_i[0, i] = - U[i] * (U[0] * Y[i, 0] * np.cos(theta[i] - theta[0])
                              + U[1] * Y[i, 1] * np.cos(theta[i] - theta[1])
                              + U[2] * Y[i, 2] * np.cos(theta[i] - theta[2]))

    # print('P_i=', P_i)
    # print('Q_i=', Q_i)
    DeltaP = P - P_i[0, 0:2]
    DeltaQ = Q - Q_i[0, 0]
    # print(DeltaP, DeltaQ)
    P = P_i[0, 0:2]
    Q = Q_i[0, 0]

    if (DeltaP[0] < 0.0001) & (DeltaP[1] < 0.0001) & (DeltaQ < 0.0001):
        T = 0

    # 形成J矩阵
    H = np.zeros([2, 2])
    N = np.zeros([2, 1])
    M = np.zeros([1, 2])
    L = np.zeros([1, 1])

    H[0, 0] = U[0] * U[0] * Y[0, 0] + Q_i[0, 0]
    H[1, 1] = U[1] * U[1] * Y[1, 1] + Q_i[0, 1]
    H[0, 1] = U[0] * U[1] * Y[0, 1] * np.cos(theta[0] - theta[1])
    H[1, 0] = U[1] * U[0] * Y[1, 0] * np.cos(theta[1] - theta[0])

    N[0, 0] = - P_i[0, 0]
    N[1, 0] = - U[1] * U[0] * Y[1, 0] * np.sin(theta[1] - theta[0])

    M[0, 0] = - P_i[0, 0]
    M[0, 1] = U[0] * U[1] * Y[0, 1] * np.sin(theta[0] - theta[1])

    L[0, 0] = U[0] * U[0] * Y[0, 0] - Q_i[0, 0]
    # print(H, N)
    temp1 = np.concatenate((H, N), axis=1)
    # print(temp1)
    temp2 = np.concatenate((M, L), axis=1)
    # print(temp2)
    J = np.concatenate((temp1, temp2), axis=0)

    # 计算修正量
    y = np.append(DeltaP, DeltaQ)
    x = - np.linalg.inv(J).dot(y)
    # print(x, y)

    # 修正变量
    theta[0] = theta[0] + x[0]
    theta[1] = theta[1] + x[1]
    U[0] = U[0] + U[0] * x[2]

    # print(J)

    iteration += 1

    # plt.scatter(iteration, U[0])

end_time = time.time()

# plt.show()
print('total_time=', end_time-start_time)
print('iteration=', iteration)
print('U=', U)
print('theta=', theta)
