import numpy as np
import time

start_time = time.time()
# 形成导纳矩阵
g1 = -10
g2 = 0.01
Y = np.array([[2*(g1+g2), -g1, -g1],
              [-g1, 2*(g1+g2), -g1],
              [-g1, -g1, 2*(g1+g2)]])


# 置初值
U = np.array([1.0, 1.05, 1.0])
theta = np.array([0.0, 0.0, 0.0])

P = np.array([-2.863, -0.6661])
Q = -1.2244

# 形成B'与B''
B1 = np.array([[20, -10],
               [-10, 20]])
B2 = 19.98

x2 = 1
x1 = np.array([1, 1])
iterations = 0
detaP = np.zeros([2, 1])
KP = KQ = 1

while (KQ == 1) | (KP == 1):

    if KP == 1:
        tempP1 = U[0] * (U[0] * Y[0, 0] * np.sin(theta[0] - theta[0])
                         + U[1] * Y[0, 1] * np.sin(theta[0] - theta[1])
                         + U[2] * Y[0, 2] * np.sin(theta[0] - theta[2]))
        tempP2 = U[1] * (U[0] * Y[1, 0] * np.sin(theta[1] - theta[0])
                         + U[1] * Y[1, 1] * np.sin(theta[1] - theta[1])
                         + U[2] * Y[1, 2] * np.sin(theta[1] - theta[2]))
        detaP[0] = P[0] - tempP1
        detaP[1] = P[1] - tempP2

        P[0] = tempP1
        P[1] = tempP2

        tempP1 = tempP1 / U[0]
        tempP2 = tempP2 / U[1]

        if (abs(tempP1) < 0.0001) & (abs(tempP2) < 0.0001):
            KP = 0

        y = np.array([tempP1, tempP2])
        x1 = np.linalg.inv(B1).dot(y)
        theta[0] += x1[0]
        theta[1] += x1[1]

    while KQ == 1:
        tempQ = - U[0] * (U[0] * Y[0, 0] * np.cos(theta[0] - theta[0])
                          + U[1] * Y[0, 1] * np.cos(theta[0] - theta[1])
                          + U[2] * Y[0, 2] * np.cos(theta[0] - theta[2]))
        detaQ = Q - tempQ
        Q = tempQ
        tempQ = detaQ / U[0]
        if tempQ < 0.0001:
            KQ = 0

        x2 = tempQ / B2

        U[0] += x2

    iterations += 1


end_time = time.time()
total_time = end_time - start_time
print('total_time=', total_time)

print('theta=', theta)
print('U=', U)
print('iterations=', iterations)


