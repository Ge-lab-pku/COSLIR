function [A, out] = COSLIR(Sigma_1_input, Sigma_2_input, mu_1_input, mu_2_input, opts)
%%  Input:
%
%   Sigma_1, Sigma_2:   two matrixs, covariance matrix;
%   mu_1, mu_2:         two vectors, mean value of the sample in two stages;
%   opts.epsilon:       a positive real number, termination condition;
%   opts.lambda:      a positive real number, controling the sparsity of A;
%   opts.eta:      a positive real number, controling the size of intercept;
%   opts.Sign:          matrix with -1, 0, 1, preknowledge of the entries in A;
%   opts.rho:           a positive real number, start value of the penalty(rho);
%   opts.silent:        0 or 1, if the information is echoed to the screen;
%   opts.start:         start value of A
%   opts.rho_first_multiply:    a positive real number, the multiplier when rho was first changed;
%   opts.rho_max:       a positive real number, the uperbound of rho;
%   opts.rho_min:       a positive real number, the lowerbound of rho;
%   opts.rho_amplify:   a positive real number no less than 1, amplify
%   coefficient of rho;
%   opts.rho_shrink:    a positive real number no less than 1, shrink
%   coefficient of rho;
%   opts.rho_output_num:    a positive integer, every iter times when
%   echcoed in screen;
%   opts.rho_update_num:    a positive integer, every iter times when
%   adjusting rho


%   opts.epsilon_2      a positive real value, threshold of diff to update
%   eta;
%   opts.update_2       a positive real value, the amplify coefficient of
%   eta
    
    Sigma_1 = Sigma_1_input/norm(Sigma_1_input - Sigma_2_input) ;
    Sigma_2 = Sigma_2_input/norm(Sigma_1_input - Sigma_2_input) ;
    mu_1 = mu_1_input/norm(mu_1_input - mu_2_input);
    mu_2 = mu_2_input/norm(mu_1_input - mu_2_input);
    
    [p, ~] = size(Sigma_1);
    
    if ~isfield(opts, 'epsilon')
        epsilon = 1e-7;
    else
        epsilon = opts.epsilon;
    end
    
    if ~isfield(opts, 'rho')
        rho = 1;
    else
        rho = opts.rho;
    end
    
    if ~isfield(opts, 'Sign')
        Sign = nan;
    else
        Sign = opts.Sign;
    end
    
    if ~isfield(opts, 'silent')
        silent = 1;
    else
        silent = opts.silent;
    end
    
    if ~isfield(opts, 'rho_first_multiply')
        rho_first_multiply = 10;
    else
        rho_first_multiply = opts.rho_first_multiply;
    end
    
    if ~isfield(opts, 'rho_max')
        rho_max = 1e4;
    else
        rho_max = opts.rho_max;
    end
    
    if ~isfield(opts, 'rho_min')
        rho_min = 1e-4;
    else
        rho_min = opts.rho_min;
    end
    
    if ~isfield(opts, 'rho_amplify')
        rho_amplify = 1.01;
    else
        rho_amplify = opts.rho_amplify;
    end
    
    if ~isfield(opts, 'rho_shrink')
        rho_shrink = 1.01;
    else
        rho_shrink = opts.rho_shrink;
    end
    
    if ~isfield(opts, 'rho_output_num')
        rho_output_num = 500;
    else
        rho_output_num = opts.rho_output_num;
    end
    
    if ~isfield(opts, 'rho_update_num')
        rho_update_num = 20;
    else
        rho_update_num = opts.rho_update_num;
    end
    
    if ~isfield(opts, 'iter_max')
        iter_max = 15000;
    else
        iter_max = opts.iter_max;
    end
    
    if ~isfield(opts, 'epsilon_2')
        epsilon_2 = 5e-3;
    else
        epsilon_2 = opts.epsilon_2;
    end
    
    if ~isfield(opts, 'update_2')
        update_2 = 1.2;
    else
        update_2 = opts.update_2;
    end
    
    if ~isfield(opts, 'eta')
        eta = 1e-4;
        eta_max = 5;
    else
        eta = min(1e-4, opts.eta);
        eta_max = opts.eta;
    end
    
    if ~isfield(opts, 'start')
        start = zeros(p);
    else
        start = opts.start;
    end
    

    lambda = opts.lambda;
    
    if ~isfield(opts, 'I')
        I = eye(p);
    else
        I = opts.I;
    end
    
    I0 = eye(p);

    B2 = start + I;
    %B1 = B2;
    A = start;
    Pi_1 = 0;
    Pi_2 = 0;
    k = 0;
    MU_11 = mu_1 * mu_1';
    MU_21 = mu_2 * mu_1';
    
    out.converge = 0;
    fval = @(B1, B2, A, Pi_1, Pi_2, mu_1, mu_2) norm(Sigma_2 - B1 * Sigma_1 * B2', 'fro')^2 + ...
           lambda * sum(sum(abs(A))) + eta*norm(mu_2 - (B1+B2)*mu_1/2)^2+ ...
           trace(Pi_1' * (B1-B2)) + rho/2 * norm(B1-B2, 'fro')^2 + ...
           trace(Pi_2' * (A+I - (B1+B2)/2)) + rho/2 * norm(A+I - (B1+B2)/2, 'fro')^2;
    P1 = Sigma_1 * (B2' * B2) * Sigma_1 + 5/8 *rho*I0 + eta/4 *MU_11;
    Q1 = Sigma_2 * B2 * Sigma_1 - Pi_1/2 + Pi_2/4 +3*rho/8 * B2 + rho/4 * (A + I) +...
       eta/2 * MU_21 - eta/4 * B2 * MU_11;

    while k < iter_max
      %fprintf('iter: %d -- 1: %5.2f; ',k, fval(B1, B2, A, Pi_1, Pi_2, mu_1, mu_2));

      
       % update B1
       B1 = Q1 / P1;
       
       
       % update B2

       P2 = Sigma_1 * (B1' * B1) * Sigma_1 + 5/8 *rho*I0 + eta/4 *MU_11;
       Q2 = Sigma_2 * B1 * Sigma_1 + Pi_1/2 + Pi_2/4 +3*rho/8 * B1 + rho/4 * (A + I) +...
           eta/2 * MU_21 - eta/4 * B1 * MU_11;

       B2 = Q2 / P2;
       
      %fprintf('3: %5.2f; ',fval(B1, B2, A, Pi_1, Pi_2, mu_1, mu_2));
       % update A
       A_temp = max(abs((B1+B2)/2 - I - Pi_2 /rho)- lambda/rho, 0) .* sign((B1+B2)/2 - I - Pi_2 /rho);
       if isnan(Sign)
           A = A_temp;
       else
           A = A_temp;
           A(A_temp .* Sign <= 0) = 0;
       end
        
       
       %update langrange multiplier
       Pi_1 = Pi_1 + rho * (B1-B2);
       Pi_2 = Pi_2 + rho * (A + I - (B1+B2)/2);
       k = k+1;
       
       P1 = Sigma_1 * (B2' * B2) * Sigma_1 + 5/8 *rho*I0 + eta/4 *MU_11;
       Q1 = Sigma_2 * B2 * Sigma_1 - Pi_1/2 + Pi_2/4 +3*rho/8 * B2 + rho/4 * (A + I) +...
       eta/2 * MU_21 - eta/4 * B2 * MU_11;
     
       AB_err = norm((B1+B2)/2 - I - A, 'fro');
       G = norm(B1*P1 - Q1, 'fro');
       BB_err = norm(B1 - B2, 'fro');

        if mod(k, rho_output_num) == 0 
            if silent == 1

               mu_err = norm(mu_2 - (B1+B2)*mu_1/2);
               fprintf('iter: %d, rho: %.8f, AB_err: %.8f, ||G||: %.8f, ||B1-B2||: %.8f, mu_err: %.8f, fval: %.8f\n',k , rho, AB_err, G, BB_err, mu_err, fval(B1, B2, A, Pi_1, Pi_2, mu_1, mu_2));
            end
        end
        
        % update eta
        if max(AB_err, G) < epsilon_2 && eta < eta_max
            eta = min(eta * update_2, eta_max);
        end
        
        if k == 20
            rho = rho * rho_first_multiply; % 第一次增大倍数需要根据lambda的增大而逐渐增大
            diff_min = 1;
        elseif mod(k, rho_update_num) == 0 && k > 100
            if AB_err < epsilon && G < epsilon && eta == eta_max
                out.converge = 1;
                break
            end
            
            if AB_err > diff_min*1.5 && AB_err > epsilon *10
                rho = min(rho * rho_amplify, rho_max);
            else
                diff_min = AB_err;
                if AB_err/G < 0.5
                    if rho > 10 * rho_min
                        rho = max(rho/(rho_shrink*1.3), rho_min);
                    else
                        rho = max(rho/rho_shrink, rho_min);
                    end
                end
            end
        end   
        
    end

    out.iter = k;
 
    out.ErrS = norm(Sigma_2 - (A+I)*Sigma_1*(A'+I), 'fro')/norm(Sigma_1 -Sigma_2, 'fro');
    out.ErrMu = norm(mu_2- I * mu_1-A*mu_1, 2)/norm(mu_2-mu_1, 2);
    out.Sparsity = length(find(abs(A)> 0))/p^2;
    out.fval = (norm(Sigma_2 - (A+I)*Sigma_1*(A'+I), 'fro')^2 + lambda * sum(sum(abs(A))) + eta * norm(mu_2 - mu_1 - A*mu_1)^2);
    fprintf('Done! Iter_num = %d, fval = %5f.\n', k, out.fval)
    fprintf('Fval: %f, ErrS: %f, ErrMu: %f, Sparsity: %f\n', out.fval, out.ErrS, out.ErrMu ,out.Sparsity);
