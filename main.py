import numpy as np
from math import sin, cos, exp, log


def mu1(t): #левая
    return 0 # третий
    #return sin(t) #второй
    #return t #первый


def mu2(t): #правая
    return 2*t+log(2) #третий
    #return cos(t) #второй
    #return 2*t #первый


def U0(x):
    return log(x+1) #третий
    #return x #второй
    #return x #первый


def f(x, t):
    return 0
    #return 0


def transpose(array):
    return np.transpose(array)


def progonka(a, b, f, h, mu1, mu2):
    t = np.arange(a, b + teta, teta)
    x = np.arange(a, b + h, h)
    nt = len(t)
    nx = len(x)
    u = np.zeros((len(x), len(t)))
    for j in range(nx):
        u[j][0] = U0(x[j])

    # for j in range(nt):
    #     u[0, j] = mu1(t[j])
    #     u[-1, j] = mu2(t[j])

    for k in range(1, nt):
        A = np.zeros(nx)
        B = np.zeros(nx)
        C = np.zeros(nx)
        D = np.zeros(nx)

        A[0] = 0
        B[0] = -1 - 2 * teta/h**2
        C[0] = teta/h**2
        D[0] = -(u[0][k-1]+teta/h**2*f(x[0], t[k]))-teta/h**2*mu1(t[k])
        for i in range(1, nx - 1):
            A[i] = teta/h**2
            B[i] = -1 - 2 * teta/h**2
            C[i] = teta/h**2
            D[i] = -(u[i][k-1]+teta/h**2*f(x[i], t[k]))
        A[nx - 1] = teta/h**2
        B[nx - 1] = -1 - 2 * teta/h**2
        C[nx - 1] = 0
        D[nx - 1] = -(u[nx-1][k-1]+teta/h**2*f(x[nx-1], t[k]))-teta/h**2*mu2(t[k])

        omega = np.zeros(nx)
        alpha = np.zeros(nx)
        beta = np.zeros(nx)

        omega[0] = B[0]
        alpha[0] = -C[0] / omega[0]
        beta[0] = D[0] / omega[0]

        for i in range(1, nx - 1):
            omega[i] = B[i] + A[i] * alpha[i - 1]
            alpha[i] = -C[i] / omega[i]
            beta[i] = (D[i] - A[i] * beta[i - 1]) / omega[i]

        omega[nx-1] = B[nx-1] + A[nx-1] * alpha[nx-2]
        alpha[nx-1] = 0
        beta[nx-1] = (D[nx-1] - A[nx-1] * beta[nx-2]) / omega[nx-1]

        u[nx-1][k] = beta[nx-1]

        for i in range(nx - 2, -1, -1):
            u[i, k] = alpha[i] * u[i + 1, k] + beta[i]

    return x, t, u


a_sb, b_sb = 0, 1
alpha1 = 1
h = 0.1
teta = h ** 2 / (2 * alpha1)
x, t, u = progonka(a_sb, b_sb, f, h, mu1, mu2)

ures = transpose(u)
with open('res.txt', 'w') as f:
    sres = ['']*len(ures)
    tres, xres = [0]*len(t), [0]*len(x)
    s = 'x; '
    for i in range(len(x)):
        s = s + str(x[i]).replace('.', ',') + ';'
    f.write(s + '\n')
    for i in range(len(ures)):
        tres[i] = str(t[i]).replace('.', ',')
        for j in range(len(ures[0])):
            sres[i] = sres[i] + str(ures[i][j]).replace('.', ',') + ';'
        f.write(tres[i] + '; ' + sres[i]+'\n')
