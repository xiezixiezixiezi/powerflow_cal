import numpy as np
import time

time_start = time.time()
g1 = -10
g2 = 0.01
Y = np.array([[2*(g1+g2), -g1, -g1],
              [-g1, 2*(g1+g2), -g1],
              [-g1, -g1, 2*(g1+g2)]])

# 置初值
e = np.array([1.0, 1.05, 1.0])
f = np.array([0.0, 0.0, 0.0])
U2 = 1.05

P = np.array([-2.8653, -0.6661])
Q = -1.2244
x = [1, 1, 1]
iteration = 0
while min(np.abs(x)) > 0.000001:
    # 计算潮流方程
    Pi, Qi = np.zeros([1, 3]), np.zeros([1, 3])
    for i in range(3):
        Pi[0, i] = e[i] * (-1) * (Y[i, 0] * f[0] + Y[i, 1] * f[1] + Y[i, 2] * f[2]) + f[i] * (Y[i, 0] * e[0] + Y[i, 1] * e[1] + Y[i, 2] * e[2])
        Qi[0, i] = f[i] * (-1) * (Y[i, 0] * f[0] + Y[i, 1] * f[1] + Y[i, 2] * f[2]) - e[i] * (Y[i, 0] * e[0] + Y[i, 1] * e[1] + Y[i, 2] * e[2])

    Usqr = e[1] * e[1] + f[1] * f[1]

    deltaP = P - Pi[0, 0:2]
    deltaQ = Q - Qi[0, 0]
    deltaUsqr = U2 * U2 - Usqr
    # print(deltaP, deltaQ, deltaUsqr)

    y = np.append(deltaP, (deltaQ, deltaUsqr), axis=0)
    # print(y)

    H, N = np.zeros([2, 2]), np.zeros([2, 2])
    M, L = np.zeros([1, 2]), np.zeros([1, 2])

    H[0, 0] = Y[0, 0] * f[0] + Y[0, 1] * f[1] + Y[0, 2] * f[2] - Y[0, 0] * f[0]
    H[1, 1] = Y[1, 0] * f[0] + Y[1, 1] * f[1] + Y[1, 2] * f[2] - Y[1, 1] * f[1]
    H[0, 1] = - Y[0, 1] * f[0]
    H[1, 0] = - Y[1, 0] * f[1]

    N[0, 0] = - (Y[0, 0] * e[0] + Y[0, 1] * e[1] + Y[0, 2] * e[2]) + Y[0, 0] * e[0]
    N[1, 1] = - (Y[1, 0] * e[0] + Y[1, 1] * e[1] + Y[1, 2] * e[2]) + Y[1, 1] * e[1]
    N[0, 1] = Y[0, 1] * e[0]
    N[1, 0] = Y[1, 0] * e[1]

    M[0, 0] = - (Y[0, 0] * e[0] + Y[0, 1] * e[1] + Y[0, 2] * e[2]) + Y[0, 0] * e[0]
    M[0, 1] = Y[0, 1] * e[0]

    L[0, 0] = Y[0, 0] * f[0] + Y[0, 1] * f[1] + Y[0, 2] * f[2] + Y[0, 0] * f[0]
    L[0, 1] = Y[0, 1] * f[0]

    R = np.array([[0, -2*e[0]]])

    S = np.array([[0, -2*f[0]]])

    temp1 = np.concatenate((H, N), axis=1)
    temp2 = np.concatenate((M, L), axis=1)
    temp3 = np.concatenate((R, S), axis=1)
    J = np.concatenate((temp1, temp2, temp3), axis=0)
    print(J)

    if iteration == 0:
        J1 = J

    x = - np.linalg.inv(J).dot(y)

    e[0] += x[0]
    e[1] += x[1]
    f[0] += x[2]
    f[1] += x[3]
    iteration += 1

time_end = time.time()
time_sum = time_end - time_start
print('time=', time_sum)
print('iterations=', iteration)
print('e=', e)
print('f=', f)

# print(J1)
