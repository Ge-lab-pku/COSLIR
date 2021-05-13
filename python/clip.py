## Tune the clip threshold and the confidence using EPR
from COSLIR import *
import pandas as pd


res = np.load('data/human/hSTR_656.npy')


threshold_list = [0, 1e-3, 3e-3, 6e-3, 1e-2]
conf_list = [0.5, 0.6, 0.7]



bootstrap_num, dim, _ = res.shape

print('bootstrap, dim', bootstrap_num, dim)

# load data
X = np.genfromtxt('data/hSTRING/ExpressionData5.csv', delimiter=',')
X = X[1:,1:]
Y = np.genfromtxt('data/hSTRING/ExpressionData6.csv', delimiter=',')
Y = Y[1:,1:]

X = X.T
Y = Y.T

print('sample size: X, ', X.shape, 'Y, ', Y.shape)

Expre = pd.read_csv('data/hSTRING/ExpressionData1.csv')
# Expre = np.load('data/mouse/mSTR_1-2_Lambda6.npy')
name = Expre['X'].str.upper()
name = pd.DataFrame(name)



# rescale the final estimator
def rescale(A, X, Y):
    mu1 = np.mean(X, axis=0)
    mu2 = np.mean(Y, axis=0)
    diff_mu = np.abs(mu2 - mu1)
    diff_mu = diff_mu[:, np.newaxis]
    mu1 = mu1[:, np.newaxis]
    Coef = np.matmul(diff_mu, mu1.T)
    return np.abs(A) * Coef


# Change the coefficient matrix A into the path table format
def MatToPath(A):
    pos = np.nonzero(A)
    source_ind = pos[1]
    target_ind = pos[0]
    # num = source_ind.shape[0]
    source = pd.DataFrame(name.iloc[source_ind])
    target = pd.DataFrame(name.iloc[target_ind])
    value = pd.DataFrame(A[target_ind, source_ind])
    source.index = value.index
    target.index = value.index
    Path = pd.concat([source, target, value], axis=1)
    Path.columns = ['Gene1', 'Gene2', 'EdgeWeight']
    return Path




ind1 = 0
ind2 = 0

for threshold in threshold_list:
    ind1 += 1
    A_sum = np.zeros((dim, dim))
    A_tot = np.zeros((dim, dim))
    for i in range(bootstrap_num):
        # print('iteration: {:d}'.format(i))
        A = res[i, :, :]
        A_sum[A > threshold] += 1
        A_sum[A < -threshold] -= 1
        A[np.abs(A) < threshold] = 0
        A_tot += A
    A_sum = A_sum / bootstrap_num
    # np.save('Aconf_thre_7_'+str(ind1), A_sum)
    A_tot = A_tot / bootstrap_num
    for conf in conf_list:
        ind2 += 1
        print("threshold: {}, confidence: {}".format(threshold, conf))
        A_tot[np.abs(A_sum) < conf] = 0
        A_final = rescale(A_tot, X, Y)
        A_table = MatToPath(A_final)
        print(A_table.shape)
        outFile = open('rankedEdges_thre' + str(ind1) + '_conf' + str(ind2) + '.csv', 'w')
        outFile.write('Gene1' + '\t' + 'Gene2' + '\t' + 'EdgeWeight' + '\n')

        count = 0
        for idx, row in A_table.sort_values('EdgeWeight', ascending=False).iterrows():
            count += 1
            outFile.write('\t'.join([row['Gene1'], row['Gene2'], str(row['EdgeWeight'])]) + '\n')
            if count % 100 == 0:
                print('count: ', count)

        outFile.close()
