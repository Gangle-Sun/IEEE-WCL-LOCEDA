%* ---------------------------------------------------------------------------------------------------------------------
%*   Created by Gangle Sun on 19 October, 2024
%*   email: sungangleseu@gmail.com / sungangle@seu.edu.cn
%*   Institute: Southeast University, China.
%* ---------------------------------------------------------------------------------------------------------------------
%*   This code implements the LOw-Coherence sEquence Design Algorithm (LOCEDA) from the following publication [1], which
%*   can generate sequences with (i) low-coherence, (ii) arbitrary lengths, (iii) any number of sequences, (iv) support
%*   for adaptable subcarrier assignments in orthogonal frequency-division multiple access (OFDMA) systems, and (v)
%*   compliance with user-defined PAPR constraints.
%*   [1] G. Sun, W. Wang, W. Xu, and C. Studer, "Low-Coherence Sequence Design Under PAPR Constraints," IEEE Wireless
%*   Commun. Lett., 2024.
%* ---------------------------------------------------------------------------------------------------------------------
%*   This paper is also available on https://arxiv.org/abs/2407.21400
%*   If you find our code and paper helpful, please cite our work. Thank you very much! ^o^
%* ---------------------------------------------------------------------------------------------------------------------
%%

function [P_best,mu_min,mu_P_per_iter, mu_min_per_iter] = LOCEDA(L, N, N_S, c, Gamma_PAPR, K, K_1, gamma, rho, I_MAX, tau_seq, tau_PAPR)

%% Initialization
% Show the status of optimization
clc; fprintf('Initialization progress: %.2f%%\n', 0);   

% Initialize algorithmic parameters
R_bound = sqrt(0.5*(1-sqrt((N-L)/(L*(N-1)))));
R_seq = R_bound;
mu_min = 1;
mu_min_previous = 1;

% We stop LOCEDA either after I_MAX iterations or if mu_P_min has not been updated for more than counter_stop_max consecutive iterations.
counter_stop = 0;        % Counter to stop the LOCEDA
counter_stop_max = 500;

% Initialize the sequence matrix P: L x N with each column being the sequence p_n and calculate the current mu_P
P = complex(randn(L,N),randn(L,N));
P = P./sqrt(sum(abs(P).^2));
G = abs(P'*P);
mu_P = max(G(~logical(eye(size(G)))));

% Initialize the PAPR matrix W: L x N_S with each column being the PAPR vector w_s
W = 1/sqrt(L)*exp(-2*1i*pi*(c-1).'*(0:(N_S-1))/N_S);

% Initialize matrices to record data during LOCEDA
mu_P_per_iter = zeros(I_MAX,1);
mu_min_per_iter = zeros(I_MAX,1);

%% Run
for iter = 1:I_MAX
    %% [LINE 4 OF ALG.1] Update sequence matrix P for K rounds
    % "L <= Gamma_PAPR" means there is no PAPR constraint
    if (iter > 1) || (L <= Gamma_PAPR)        
        for i_K = 1:K
            % Calculate sequence-sequence position matrix D: (N-1) x N x L, which is composed of all sequence-sequence position vectors d_{n,m}
            P_tmp = P'*P; % P_tmp: N x N
            D = permute(repmat(P,1,1,N), [3 2 1]) - permute(repmat(P,1,1,N),[2 3 1]) .* repmat(P_tmp./abs(P_tmp),1,1,L);
            D = reshape(D(repmat(~eye(N), [1, 1, L])), [N-1, N, L]);

            % Calculate D_norm: (N-1) x N to find sequence-sequence hypersphere collisions
            D_norm = sqrt(sum(abs(D).^2,3));
            mask_D_norm = D_norm < 2*R_seq; % mask_D_norm: (N-1) x N

            % Calculate U: N x L, which is composed of all displacement vectors u_n
            U = squeeze(sum(repmat((2*R_seq-D_norm).*mask_D_norm./D_norm,1,1,L).*D,1)).';

            % Update P to reduce their coherence
            P = P + tau_seq * U;
            P = NormalizeMatrix(P);
        end
    end

    %% [LINE 5 OF ALG.1] Update sequence matrix P to satisfy the PAPR constraints
    % "L > Gamma_PAPR" means there are PAPR constraints
    if L > Gamma_PAPR        
        while 1
            % Calculate sequence-PAPR position matrix bar_D: N_S x N x L, which is composed of all sequence-PAPR position vectors bar_d_{n,s}
            bar_D = permute(repmat(P,1,1,N_S), [3 2 1]) - permute(repmat(W,1,1,N),[2 3 1]) .*repmat(W'*P./abs(W'*P),1,1,L);

            % Calculate bar_D_norm: N_S x N to find sequence-PAPR hypersphere collisions
            bar_D_norm = sqrt(sum(abs(bar_D).^2,3));
            mask_bar_D_norm = bar_D_norm < sqrt(2*(1-sqrt(Gamma_PAPR/L))); % mask_bar_D_norm: N_S x N

            % Calculate bar_U: N x L, which is composed of all displacement vectors bar_u_n
            bar_U = squeeze(sum(repmat((sqrt(2*(1-sqrt(Gamma_PAPR/L)))-bar_D_norm).*mask_bar_D_norm./bar_D_norm,1,1,L).*bar_D,1)).';

            % Exit when all sequences satisfy the PAPR constraints
            if sum(mask_bar_D_norm,'all') == 0
                if iter == 1
                    clc
                    fprintf('Initialization completed')
                end
                break
            end

            % Update P to satisfy PAPR constraints
            P = P + tau_PAPR * NormalizeMatrix(bar_U);
            P = NormalizeMatrix(P);

            % Show the status of LOCEDA
            if iter == 1
                clc
                fprintf('Initialization progress: %.2f%%\n', 100*(sum(sum(mask_bar_D_norm,1)==0)/N))
            else
                clc;
                fprintf(['Length L = ' num2str(L) ', ' 'Number N = ' num2str(N) ', ' 'PAPR threshold Gamma_PAPR = ' num2str(Gamma_PAPR) '.\n'])
                fprintf(['step sizes: tau_seq = ' num2str(tau_seq) ', tau_PAPR = ' num2str(tau_PAPR) ', gamma = ' num2str(gamma) '.\n'])
                fprintf(['Iteration: ' num2str(iter) '/' num2str(I_MAX) '.\n'])
                fprintf(['Coherence: ' num2str(mu_P) '.\n'])
                fprintf(['Difference to the Welch lower bound: ' num2str((mu_P-sqrt((N-L)/(L*(N-1))))) '.\n'])
                fprintf('Handling the PAPR constraints: %.2f%%\n', 100*(sum(sum(mask_bar_D_norm,1)==0)/N))
            end
        end
    end
    G = abs(P'*P);
    mu_P = max(G(~logical(eye(size(G)))));

    % Print the status of LOCEDA
    clc;
    fprintf(['Length L = ' num2str(L) ', ' 'Number N = ' num2str(N) ', ' 'PAPR threshold Gamma_PAPR = ' num2str(Gamma_PAPR) '.\n'])
    fprintf(['step sizes: tau_seq = ' num2str(tau_seq) ', tau_PAPR = ' num2str(tau_PAPR) ', gamma = ' num2str(gamma) '.\n'])
    fprintf(['Iteration: ' num2str(iter) '/' num2str(I_MAX) '.\n'])
    fprintf(['Coherence: ' num2str(mu_P) '.\n'])
    fprintf(['Difference to the Welch lower bound: ' num2str((mu_P-sqrt((N-L)/(L*(N-1))))) '.\n'])

    %% [LINE 6 OF ALG.1] Update R_seq and Gamma_mut
    Gamma_mut = 1-2*R_seq^2;
    if mu_P < mu_min                                 % Case 1: mu_P < mu_min
        R_seq = min(R_seq + gamma, R_bound);
        counter_stop = 0;
    else
        if mu_P <= Gamma_mut                         % Case 2: mu_min <= mu_P <= Gamma_mut
            R_seq = R_bound;
        else                                         % Case 3: mu_P > Gamma_mut
            R_seq = R_seq - gamma;
        end
    end
    % Update the counter_stop
    counter_stop = counter_stop + 1;
    if counter_stop > counter_stop_max
        break
    end

    %% [LINE 7 OF ALG.1] Update mu_min and P_best
    if mu_P < mu_min
        mu_min = mu_P;
        P_best = P;
    end
    mu_P_per_iter(iter) = mu_P;
    mu_min_per_iter(iter) = mu_min;

    %% [LINE 8 OF ALG.1] Update step sizes tau_seq and tau_PAPR every K_1 iterations
    if mod(iter, K_1) == 0
        if mu_min_per_iter(iter-1) < mu_min_previous
            tau_seq = (1 + rho) * tau_seq;
            tau_PAPR = (1 + rho) * tau_PAPR;
        else
            tau_seq = (1 - rho) * tau_seq;
            tau_PAPR = (1 - rho) * tau_PAPR;
        end
        mu_min_previous = mu_min_per_iter(iter-1);
    end

end

end