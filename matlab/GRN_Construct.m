function GRN_Construct(input1, input2, Lambda, Repeat, outfile, seed, opts)

%opts = struct();
%opts.epsilon = 1e-7;
opts.rho = 1;
%opts.Sign = nan;
%opts.silent = 0;
opts.rho_first_multiply = 10;
opts.rho_max = 1e4;
opts.rho_min = 1e-4;
opts.rho_amplify = 1.01;
opts.rho_shrink = 1.02;
%opts.rho_output_num = 500;
%opts.rho_update_num = 20;
%opts.iter_max = 15000;
opts.epsilon_2 = 5e-3;
opts.update_2 = 1.2;
%opts.eta = 5;
dim1 = size(input1);
dim2 = size(input2);
N1 = dim1(2);
N2 = dim2(2);

rng(seed);
Rec_2 = cell(Repeat ,1);
for i = 1:Repeat
    opts_temp = opts;
    r1 = randi(N1, N1, 1);
    r2 = randi(N2, N2, 1);
    Sample1_temp = input1(:,r1);
    Sample2_temp = input2(:,r2);
    Cov1 = cov(Sample1_temp');
    Cov2 = cov(Sample2_temp');
    opts.start = zeros(dim1(1));
    Mean1 = mean(Sample1_temp')';
    Mean2 = mean(Sample2_temp')';
    Rec_2{i} = cell(length(Lambda), 1);
    for j = 1:length(Lambda)
        opts_temp.lambda = Lambda(j);
        [A1, out] = COSLIR(Cov1, Cov2, Mean1, Mean2, opts_temp);
        Rec_2{i}{j} = {A1, out, Mean1, Mean2};
    end
end

save(outfile, 'Rec_2');
end

