## An example script to run COSLIR
from COSLIR import *

input1 = '../data/example/ExpressionData1.csv'
input2 = '../data/example/ExpressionData2.csv'
Lambda = 1e-7
eta = 5
epsilon = 1e-4
screen = True
iter_max = 200000
bootnum = 2
checkiter = 10
iters_per_screen = 100

X = np.genfromtxt(input1, delimiter=',')
X = X[1:, 1:]
Y = np.genfromtxt(input2, delimiter=',')
Y = Y[1:, 1:]

X = X.T
Y = Y.T

print('sample size: X, ', X.shape, 'Y, ', Y.shape)
opts = {
    'eta': eta,
    'epsilon': epsilon,
    'screen': screen,
    'iter_max': iter_max,
    'iters_per_screen': iters_per_screen,
    'rho_update_num': 2000,
    'rho_shrink': 1.01
}
solver = COSLIR_bootstrap(X, Y, bootnum, Lambda, opts)
solver.bootstrap(checkiter=checkiter, seed=1)