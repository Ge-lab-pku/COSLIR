import sys
import argparse
from COSLIR import *


parser = argparse.ArgumentParser(description='Gene regulation netowrk inference using COSLIR')
######################################
## essential parameters

##parser.add_argument('--train', default=True, help='If set true, train the model')

parser.add_argument('--input1', help='The expression data for the first stage')

parser.add_argument('--input2', help='The expression data for the second stage')

parser.add_argument('--Lambda', type=float, help='hyper-parameter that controls the sparsity',
                    default=1e-6)
parser.add_argument('--eta', type=float, help='hyper-parameter that penalize the size of intercept',
                    default=5)
parser.add_argument('--epsilon', type=float, help='hyper-parameter that determines the terminate condition',
                    default=1e-5)
parser.add_argument('--screen', type=bool, help='determines whether the convergent behavior is echoed to \
                    the screen', default=True)
parser.add_argument('--iter_max', type=int, help='the maximal iteration number',
                    default=20000)
parser.add_argument('--bootnum', type=int, help='the bootstrapping number',
                    default=50)
parser.add_argument('--checkiter', type=int, help='every iter times in bootstrapping will make a checkpoint and save\
                      the data', default=10)

######################################
# parser.add_argument('--rho', type=float, help='start value of the penalty (rho) in COSLIR',
#                     default=1)
parser.add_argument('--iters_per_screen', type=int, help='every iter times when echcoed in screen',
                    default=20)
# parser.add_argument('--start_value', help='provides the initial value of the A matrix',
#                     default='default')

# parser.add_argument('--rho_first_multiply', type=float, help='the multiplier that changes rho in the process',
#                     default=10)
# parser.add_argument('--rho_max', type=float, help='an upper bound of rho',
#                     default=1e4)
# parser.add_argument('--rho_min', type=float, help='a lower bound of rho',
#                     default=1e-4)
# parser.add_argument('--rho_amplify', type=float, help='the amplify coefficient of rho, positive and no less than 1',
#                     default=1.02)
parser.add_argument('--rho_shrink', type=float, help='the shrink coefficient of rho, positive and no less than 1',
                    default=1.01)
parser.add_argument('--rho_update_num', type=int, help='every iter times when adjusting rho',
                    default=20)
# parser.add_argument('--epsilon_eta', type=float, help='threshold of diff to update eta',
#                     default=5e-3)
# parser.add_argument('--update_eta', type=float, help='the amplify coefficient of eta',
#                     default=1.2)




if __name__ == '__main__':
    args = parser.parse_args()
    X = np.genfromtxt(args.input1, delimiter=',')
    X = X[1:, 1:]
    Y = np.genfromtxt(args.input2, delimiter=',')
    Y = Y[1:, 1:]

    X = X.T
    Y = Y.T

    print('sample size: X, ', X.shape, 'Y, ', Y.shape)
    opts = {
        'eta': args.eta,
        'epsilon': args.epsilon,
        'screen': args.screen,
        'iter_max': args.iter_max,
        'iters_per_screen': args.iters_per_screen,
        'rho_update_num': args.rho_update_num,
        'rho_shrink': args.rho_shrink
    }
    solver = COSLIR_bootstrap(X, Y, args.bootnum, args.Lambda, opts)
    solver.bootstrap(checkiter=args.checkiter, seed=1)


