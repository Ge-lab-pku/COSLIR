%% Simple example that runs COSLIR on the real data sets
clc;
% Read inputs
input1_file = '../data/example/ExpressionData1.csv';
input2_file = '../data/example/ExpressionData2.csv';

input1 = csvread(input1_file, 1, 1);
input2 = csvread(input2_file, 1, 1);
% Specify parameters
outfile = 'output';
bootnum = 2;
Lambda = [1e-6];
seed = 'default';
opts = struct();
opts.eta = 5;
opts.epsilon = 1e-4;
opts.silent = 1;
opts.iter_max = 500000;
opts.rho_output_num = 1000;
opts.rho_update_num = 2000;


% Run oracle test%%
GRN_Construct(input1, input2, Lambda, bootnum, outfile, seed, opts)