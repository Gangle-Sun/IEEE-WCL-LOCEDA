%* ---------------------------------------------------------------------------------------------------------------------
%*   Created by Gangle Sun on 19 October, 2024
%*   email: sungangleseu@gmail.com / sungangle@seu.edu.cn
%*   Institute: Southeast University, China.
%* ---------------------------------------------------------------------------------------------------------------------
%*   This code implements the LOw-Coherence sEquence Design Algorithm (LOCEDA) from the following publication [1], which 
%*   can generate sequences with (i) low-coherence, (ii) arbitrary lengths, (iii) any number of sequences, (iv) support 
%*   for adaptable subcarrier assignments in orthogonal frequency-division multiple access (OFDMA) systems, and (v) 
%*   compliance with user-defined PAPR constraints.
%*   [1] G. Sun, W. Wang, W. Xu, and C. Studer, "Low-Coherence Sequence  Design Under PAPR Constraints," IEEE Wireless 
%*   Commun. Lett., 2024.
%* ---------------------------------------------------------------------------------------------------------------------
%*   This paper is also available on https://arxiv.org/abs/2407.21400 
%*   If you find our code and paper helpful, please cite our work. Thank you very much! ^o^
%* ---------------------------------------------------------------------------------------------------------------------
%%

clear all; clc;

%% parmeter setup
% scenario parameters
L_list = 36:12:108;               % List of suquence lengths
N_list = 200;                     % List of the numbers of sequences
N_C = 1024;                       % Number of total subcarriers in OFDMA systems
N_S = N_C;                        % Number of time-domain discrete sampled signals
Gamma_PAPR_list = [1.5 2 4 1e3];  % List of PAPR thresholds

% algorithm parameters
gamma = 1e-4;                     % Step size for updating R_seq
I_MAX = 1e4;                      % Maximum iteration number
tau_seq = 5e-2;                   % Step size for updating sequences to reduce coherence
tau_PAPR = 5e-2;                  % Step size for updating sequences to satisfy PAPR constraints
K = 5;                            % Round number to update sequences every time
K_1 = 20;                         % Update tau_seq and tau_PAPR every K_1 iterations
rho = 0.05;                       % parameter to control the speed for updating tau_seq and tau_PAPR

%% main
for i_L = 1:length(L_list)
    for i_N = 1:length(N_list)
        for i_Gamma_PAPR =1:length(Gamma_PAPR_list)

            L = L_list(i_L);      % Suquence length
            N = N_list(i_N);      % Number of sequences
            Gamma_PAPR = Gamma_PAPR_list(i_Gamma_PAPR);     % PAPR threshold
            c = 1:L;              % Index vector of assigned subcarriers

            % Run LOCEDA
            [P,mu_P_best,mu_P_per_iter, mu_min_per_iter] = LOCEDA(L, N, N_S, c, Gamma_PAPR, K, K_1, gamma, rho, I_MAX, tau_seq, tau_PAPR);
            filename = ['LOCEDA_with_L_' num2str(L) '_N_' num2str(N) '_Gamma_PAPR_' num2str(Gamma_PAPR) '.mat'];
            save(filename,'P','mu_P_best','mu_P_per_iter','mu_min_per_iter');
        end
    end
end