import numpy as np
from numpy.linalg import norm
from IPython import embed

class COSLIR:

    def __init__(self,
                 Lambda,
                 eta=5,
                 epsilon=2.5e-5,
                 rho=1,
                 screen=False,
                 iters_per_screen=20,
                 start_value='default',
                 rho_first_multiply=10,
                 rho_max=1e4,
                 rho_min=1e-4,
                 rho_amplify=1.01,
                 rho_shrink=1.01,
                 rho_update_num=20,
                 epsilon_eta=5e-3,
                 update_eta=1.2,
                 iter_max=20000):

        r'''epsilon:       a positive real number, termination condition;
            lambda:        a positive real number, controling the sparsity of A;
            eta:           a positive real number, controling the size of intercept;
            Sign:          matrix with -1, 0, 1, preknowledge of the entries in A;
            rho:           a positive real number, start value of the penalty(rho);
            screen:        a bool value, if the information is echoed to the screen;
            start_value:   start value of A
            rho_first_multiply:    a positive real number, the multiplier when rho was first changed;
            rho_max:       a positive real number, the uperbound of rho;
            rho_min:       a positive real number, the lowerbound of rho;
            rho_amplify:   a positive real number no less than 1, amplify coefficient of rho;
            rho_shrink:    a positive real number no less than 1, shrink coefficient of rho;
            rho_output_num:    a positive integer, every iter times when echcoed in screen;
            rho_update_num:    a positive integer, every iter times when adjusting rho
            epsilon_eta      a positive real value, threshold of diff to update eta;
            update_eta       a positive real value, the amplify coefficient of eta
        '''

        self.Lambda = Lambda

        self.eta_max = eta
        self.eta = min(1e-4, eta)

        self.epsilon = epsilon
        self.rho = rho
        self.screen = screen
        self.iters_per_screen = iters_per_screen
        self.start_value = start_value
        self.rho_first_multiply = rho_first_multiply
        self.rho_max = rho_max
        self.rho_min = rho_min
        self.rho_amplify = rho_amplify
        self.rho_shrink = rho_shrink
        self.rho_update_num = rho_update_num

        self.epsilon_eta = epsilon_eta
        self.update_eta = update_eta
        self.iter_max = iter_max


    def admm(self, S1, S2, m1, m2):
        sn = norm(S1 - S2, 'fro')
        self.Sigma_1 = S1 / sn
        self.Sigma_2 = S2 / sn
        mn = norm(m1 - m2)
        self.mu_1 = m1 / mn
        self.mu_2 = m2 / mn
        self.I = np.eye(S1.shape[0])
        if self.start_value == 'default':
            self.A = np.zeros(S1.shape)
        else:
            self.A = self.start_value

        self.B1 = self.A + self.I
        self.B2 = self.A + self.I
        self.Pi_1 = np.zeros(self.A.shape)
        self.Pi_2 = np.zeros(self.A.shape)
        self.MU_11 = np.matmul(self.mu_1.reshape(-1,1), self.mu_1.reshape(1,-1))
        self.MU_21 = np.matmul(self.mu_2.reshape(-1,1), self.mu_1.reshape(1,-1))
        self.P1 = self.Sigma_1 @ (self.B2.T @ self.B2) @ self.Sigma_1 + 5/8 * self.rho * self.I  + self.eta/4 * self.MU_11
        self.Q1 = self.Sigma_2 @ self.B2 @ self.Sigma_1 - self.Pi_1/2 + self.Pi_2/4 + 3/8 * self.rho * self.B2 + \
                  self.rho/4 * (self.A + self.I) + self.eta/2 * self.MU_21 - self.eta/4 * self.B2 @ self.MU_11

        #print(self.MU_11.sum(), self.MU_21.sum(), self.P1.sum(), self.Q1.sum())

        output = {}
        output['resume'] = False
        output['if_converge'] = False
        for k in range(self.iter_max):

            #print('1 {}'.format(self.fval()))
            self.B1 = self.Q1 @ np.linalg.inv(self.P1)
            #print('2 {}'.format(self.fval()))

            self.P2 = self.Sigma_1 @ (self.B1.T @ self.B1) @ self.Sigma_1 + 5/8 * self.rho * self.I  + self.eta/4 * self.MU_11
            self.Q2 = self.Sigma_2 @ self.B1 @ self.Sigma_1 + self.Pi_1/2 + self.Pi_2/4 + 3/8 * self.rho * self.B1 + \
                  self.rho/4 * (self.A + self.I) + self.eta/2 * self.MU_21 - self.eta/4 * self.B1 @ self.MU_11

            self.B2 = self.Q2 @ np.linalg.inv(self.P2)
            #print('3 {}'.format(self.fval()))

            self.A = (self.B1+self.B2)/2 - self.I - self.Pi_2/self.rho
            self.A[np.abs(self.A) < self.Lambda/self.rho] = 0
            self.A[self.A > 0] -= self.Lambda/self.rho
            self.A[self.A < 0] += self.Lambda/self.rho

            #print('4 {}'.format(self.fval()))

            self.Pi_1 = self.Pi_1 + self.rho * (self.B1 - self.B2)
            self.Pi_2 = self.Pi_2 + self.rho * (self.A + self.I - (self.B1 + self.B2)/2)

            #print('5 {}'.format(self.fval()))

            self.P1 = self.Sigma_1 @ (self.B2.T @ self.B2) @ self.Sigma_1 + 5/8 * self.rho * self.I  + self.eta/4 * self.MU_11
            self.Q1 = self.Sigma_2 @ self.B2 @ self.Sigma_1 - self.Pi_1/2 + self.Pi_2/4 + 3/8 * self.rho * self.B2 + \
                  self.rho/4 * (self.A + self.I) + self.eta/2 * self.MU_21 - self.eta/4 * self.B2 @ self.MU_11

            AB_err = norm((self.B1+self.B2)/2 - self.I - self.A, 'fro')
            G = norm(self.B1 @ self.P1 - self.Q1, 'fro')
            BB_err = norm(self.B1 - self.B2, 'fro')

            if (k+1) % self.iters_per_screen == 0 and self.screen:
                mu_err = norm(self.mu_2 - (self.B1+self.B2)@self.mu_1/2)
                print('iter: {:d}, rho: {:.8f}, AB_err: {:.8f}, ||G||: {:.8f}, ||B1-B2||: {:.8f}, ||mu_err||: {:.8f}, fval: {:.8f}'.format(k, self.rho, AB_err, G, BB_err, mu_err, self.fval()))

            if max(AB_err, G) < self.epsilon_eta and self.eta < self.eta_max:
                self.eta = min(self.eta * self.update_eta, self.eta_max)

            if k == 19:
                self.rho *= self.rho_first_multiply
                diff_min = 1
                # if AB_err > 0.005 or 
                # if G > 1:
                #     output['resume'] = True
                #     break
            elif (k+1) % self.rho_update_num == 0 and k > 100:
                if AB_err < self.epsilon and G < self.epsilon and self.eta == self.eta_max:
                    output['if_converge'] = True
                    break

                elif AB_err > diff_min*1.5 and AB_err > self.epsilon *10:
                    self.rho = min(self.rho * self.rho_amplify, self.rho_max)
                else:
                    diff_min = AB_err
                    if AB_err / G < 0.5:
                        if self.rho > 10 * self.rho_min:
                            self.rho = max(self.rho/(self.rho_shrink * 1.3), self.rho_min)
                        else:
                            self.rho = max(self.rho/self.rho_shrink, self.rho_min)

       
        output['Final_iter'] = k
        output['ErrS'] = norm(self.Sigma_2 - (self.A + self.I) @ self.Sigma_1 @ (self.A.T + self.I), 'fro')
        output['ErrMu'] = norm(self.mu_2 - (self.I + self.A) @ self.mu_1)
        output['Sparsity'] = (self.A != 0).sum() / self.A.shape[0]**2
        output['Fval'] = norm(self.Sigma_2 - (self.A + self.I) @ self.Sigma_1 @ (self.A.T + self.I), 'fro')**2 + \
                         self.Lambda * np.abs(self.A).sum() + self.eta_max * norm(self.mu_2 - self.mu_1 - self.A @ self.mu_1)**2
        print('Done! Iter_num= {:d}, fval = {:.5f}.'.format(output['Final_iter'], output['Fval']))
        print('Fval: {:f}, ErrS: {:f}, ErrMu: {:f}, Sparsity: {:f}'.format(output['Fval'], output['ErrS'], output['ErrMu'] ,output['Sparsity']))

        return self.A, output

    def fval(self):
        val = norm(self.Sigma_2 - self.B1 @ self.Sigma_1 @ self.B2.T, 'fro')**2 + \
            self.Lambda * np.abs(self.A).sum() + self.eta * norm(self.mu_2 - (self.B1+self.B2)/2 @ self.mu_1)**2 + \
            (self.Pi_1.T @ (self.B1 - self.B2)).trace() + self.rho/2 * norm(self.B1 - self.B2, 'fro')**2 + \
            (self.Pi_2.T @ (self.A + self.I - (self.B1 + self.B2)/2)).trace() + self.rho/2 * norm(self.A+self.I - (self.B1+self.B2)/2, 'fro')**2

        return val

    def compare(self, A0, threshold):
        A_temp = self.A.copy()
        A_temp[abs(self.A)<threshold] = 0
        error = norm(A_temp-A0, 'fro')/norm(A0, 'fro')
        precision = (A_temp * A0 > 0).sum() / (A_temp != 0).sum()
        recall = (A_temp * A0 > 0).sum() / (A0 != 0).sum()
        print('Threshold: {}, error: {}, precsion: {}, recall: {}'.format(threshold, error, precision, recall))
        return precision, recall

    def threshold(self, thres_list=[0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-5]):
        thres_list.sort()
        A = self.A.copy()
        for thres in thres_list:
            A[abs(A)<thres] = 0
            sparsity = (A != 0).sum() / A.shape[0]**2
            fval = norm(self.Sigma_2 - (A + self.I) @ self.Sigma_1 @ (A.T + self.I), 'fro')**2 + \
                        self.eta_max * norm(self.mu_2 - self.mu_1 - A @ self.mu_1)**2
            print('thres: {}, sparsity: {}, fval: {}'.format(thres, sparsity, fval))


class COSLIR_bootstrap:

    def __init__(self, X, Y, bootstrap_num, Lambda, opts):

        self.X = X
        self.Y = Y
        self.bootstrap_num = bootstrap_num
        self.solver = COSLIR(Lambda, eta=opts['eta'], epsilon=opts['epsilon'],
                             screen=opts['screen'], iters_per_screen=opts['iters_per_screen'],
                             iter_max=opts['iter_max'], rho_update_num=opts['rho_update_num'],
                             rho_shrink=opts['rho_shrink'])

    def sample(self, X, seed=1):
        np.random.seed(seed)
        index = list(np.random.choice(X.shape[0], size=X.shape[0], replace=True))
        X_temp = X[index, :]
        return X_temp

    def bootstrap(self, checkiter=10, seed=1):
        self.rec = []
        self.errS = []
        self.errmu = []
        self.spar = []
        errS_sum = 0
        errmu_sum = 0
        sparsity = 0 
        i = 0
        while i < self.bootstrap_num:
            self.solver.rho = 1
            print("itertaion: ", i)
            i += 1
            X_temp = self.sample(self.X, seed+i)
            Y_temp = self.sample(self.Y, seed+i*2)
            cov1 = np.cov(X_temp.T)
            cov2 = np.cov(Y_temp.T)
            mean1 = np.mean(X_temp, axis=0)
            mean2 = np.mean(Y_temp, axis=0)
            A, out = self.solver.admm(cov1, cov2, mean1, mean2)
            # if out['resume'] == True:
            #    i -= 1
            #    continue
            self.rec.append(A)
            self.errS.append(out['ErrS'])
            self.errmu.append(out['ErrMu'])
            self.spar.append(out['Sparsity'])
            errS_sum += out['ErrS']
            errmu_sum += out['ErrMu']
            sparsity += out['Sparsity']
            if out['if_converge'] != True:
                np.save('hCELL_356', self.rec)
                #np.save('mNon_734_errs', self.errS)
                #np.save('mNon_734_errmu', self.errmu)
                #np.save('mNon_734_spar', self.spar)
                break
            if i%checkiter == 0:
                np.save('hCELL_356', self.rec)
        score = [0,0,0]
        score[0] = errS_sum/self.bootstrap_num
        score[1] = errmu_sum/self.bootstrap_num
        score[2] = sparsity/self.bootstrap_num
        print('ErrS: {:f}, ErrMu: {:f}, Sparsity: {:f}'.format(score[0], score[1], score[2]))
        np.save('hCELL_356', self.rec)              
        #np.save('mNon_734_errs', self.errS)
        #np.save('mNon_734_errmu', self.errmu)
        #np.save('mNon_734_spar', self.spar)


    def eval_threshold(self, threshold, A):
        I = np.eye(A.shape[0])
        cov1 = np.cov(self.X.T)
        cov2 = np.cov(self.Y.T)
        diff_cov = norm(cov2-cov1)
        cov1 /= diff_cov
        cov2 /= diff_cov
        mu1 = np.mean(self.X, axis=0)
        mu2 = np.mean(self.Y, axis=0)
        diff_mu = norm(mu2-mu1)
        mu1 /= diff_mu
        mu2 /= diff_mu

        A_temp = A
        A_temp[np.abs(A_temp) < threshold] = 0

        loss =  norm((A_temp+I)@cov1@(A_temp+I).T-cov2)**2 + self.solver.Lambda * norm(mu2-mu1-A_temp@mu1)**2
        sparsity = (A_temp != 0).sum() / A.shape[0]**2

        return loss, sparsity
        

    def summary(self, threshold):
        self.A_sum = np.zeros((self.X.shape[1], self.X.shape[1]))
        for i in range(len(self.rec)):
            self.A_sum[self.rec[i] > threshold] += 1
            self.A_sum[self.rec[i] < -threshold] -= 1

        self.A_sum /= len(self.rec)

    def confidence(self, conf, A0):
        A_temp = self.A_sum.copy()
        A_temp[abs(self.A_sum) < conf] = 0
        precision = (A_temp * A0 > 0).sum() / (A_temp != 0).sum()
        recall = (A_temp * A0 > 0).sum() / (A0 != 0).sum()
        print('Confidence: {}, precsion: {}, recall: {}'.format(conf, precision, recall))
        return precision, recall

def oracle(dim, sample_size, seed):
    s = 0.1
    np.random.seed(seed)
    P = np.random.normal(size=(dim, dim)) + np.eye(dim)
    S1 = P @ np.diag(np.exp(np.arange(1/dim, 1+1/dim, 1/dim))) @ P.T

    index = list(np.random.choice(dim**2, size=int(np.floor(dim**2 * s)), replace=False))
    A0 = np.zeros(dim**2)
    A0[index] = np.random.normal(size=len(index))/10
    A0 = A0.reshape(dim, dim)
    m1 = np.random.normal(size=dim) * 100
    m2 = (np.eye(dim)+A0) @ m1 + np.ones(dim) * 0.1
    S2 = (np.eye(dim)+A0) @ S1 @ (np.eye(dim)+A0).T

    X = P @ np.diag(np.exp(np.arange(1/dim, 1+1/dim, 1/dim)/2)) @ np.random.normal(size=(dim, sample_size))
    X = (X + m1.reshape(dim,1)).T
    Y = (np.eye(dim)+A0) @ P @ np.diag(np.exp(np.arange(1/dim, 1+1/dim, 1/dim)/2)) @ np.random.normal(size=(dim, sample_size))
    Y = (Y + m2.reshape(dim,1)).T

    return S1, S2, m1, m2, A0, X, Y

def main():
    import pickle
    # from Norm_test import oracle
    S1, S2, m1, m2, A0 = oracle(100, A_norm=1, seed=1)
    solver = COSLIR(1e-6, screen=True)
    solver.admm(S1, S2, m1, m2)
    solver.threshold()
    solver.compare(A0, [0.0])


if __name__ == '__main__':
    main()
