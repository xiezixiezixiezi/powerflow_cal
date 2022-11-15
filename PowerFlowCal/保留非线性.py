import numpy as np
from matplotlib import pyplot as plt
import time

time_start = time.time()
plt.style.use('ggplot')


def cal_power_flow(e, f):
    # 计算潮流方程
    Pi, Qi = np.zeros([1, 3]), np.zeros([1, 3])
    for i in range(3):
        Pi[0, i] = e[i] * (-1) * (Y[i, 0] * f[0] + Y[i, 1] * f[1] + Y[i, 2] * f[2]) + f[i] * (
                    Y[i, 0] * e[0] + Y[i, 1] * e[1] + Y[i, 2] * e[2])
        Qi[0, i] = f[i] * (-1) * (Y[i, 0] * f[0] + Y[i, 1] * f[1] + Y[i, 2] * f[2]) - e[i] * (
                    Y[i, 0] * e[0] + Y[i, 1] * e[1] + Y[i, 2] * e[2])
    Usqr = e[1] * e[1] + f[1] * f[1]
    return Pi, Qi, Usqr


J = np.array([[0, 0, -20.5, 10],
              [0, 0, 10.5, -20],
              [-20.5, 10, 0, 0],
              [0, -2, 0, 0]])

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
x = [0, 0, 0, 0]
iteration = 0

while (abs(min(x)) > 0.00001) | (iteration == 0):

    Pi, Qi, U = cal_power_flow(e, f)
    dx1 = np.array(x)
    dx2 = np.array(x)
    np.insert(dx1, -1, 0)
    np.insert(dx2, -1, 0)
    # 在最后一位插入0
    Pi1, Qi1, U1 = cal_power_flow(dx1, dx2)

    deltaP = P - Pi[0, 0:2] + Pi1[0, 0:2]
    deltaQ = Q - Qi[0, 0] + Qi1[0, 0]
    deltaUsqr = U2 * U2 - U + U1
    # print(deltaP, deltaQ, deltaUsqr)

    y = np.append(deltaP, (deltaQ, deltaUsqr), axis=0)
    # print(y)

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

