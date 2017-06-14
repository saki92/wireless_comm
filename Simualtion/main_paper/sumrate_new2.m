clear;
clc;

Nt = 2; %no of transmit antennas
K = 2; %no of receivers
Q = 1; %no of antennas per receiver
M = 5; %no of error samples
MV = 20; %no of error samples for rate validation
sigma = 1; %noise covariance
SNR = -5:2:35; %dB
PtL = 10.^(SNR/10); %total transmit power
noofit = 10; %no of validation iteration
error_var = 0.3; %alpha factor


%%%%%%%%%%% Channel %%%%%%%%%
Hcap = ( 1/sqrt(2) ) * ( randn(Nt,K) + 1i*randn(Nt,K) )*sqrt(1- error_var); %estimate
Htil = ( 1/sqrt(2) ) * ( randn(Nt,K,M,noofit) + 1i*randn(Nt,K,M,noofit) )*sqrt(error_var); %error samples
H = Hcap + Htil; %actual

HtilVal = ( 1/sqrt(2) ) * ( randn(Nt,K,MV) + 1i*randn(Nt,K,MV) )*sqrt(error_var);
HVal = Hcap + HtilVal;
%%%%%%%%%%% Channel %%%%%%%%%

P = ( 1 / sumsqr(abs(Hcap')) * Hcap' ).'; %initial Tx vectors for all samples

AvRate = zeros(length(SNR),1);
for p = 1:length(PtL) %for each SNR values
    Pt = PtL(p);
    for ite = 1:noofit %main loop for averaging Sum Rate

        a = 1;
        n = 0;
        SVal = [];
        while a == 1 % Loop for convergance of Tx matrix
            T = zeros(M*Nt, K); G = zeros(M*Nt, K); U = zeros(M*Nt, K); v = zeros(Nt*M,K);
            t = zeros(M*Nt, K); Psi = zeros(M*Nt,K,Nt,Nt); f = zeros(M*Nt,K,Nt);
            for k = 1:K
                for m = 1:M
                    Ik = 0;
                    for i = 1:K
                        if i ~= k
                            Ik = Ik + abs( H(:,k,m,ite)' * P(:,i) )^2 + sigma; %Interfearence power
                        end
                    end
                    T(m,k) = abs( H(:,k,m,ite)' * P(:,k) )^2 + Ik; %Average power
                    G(m,k) = P(:,k)' *H(:,k,m,ite) / T(m,k); %Equalizer matrix
%                     U(m,k) = inv( T(m,k) \ Ik); %MMSE weights
                    U(m,k) = T(m,k) / Ik; %MMSE weights
                    t(m,k) = U(m,k) * abs( G(m,k) )^2;
                    Psi(m,k,:,:) = t(m,k) * H(:,k,m,ite) * H(:,k,m,ite)';
                    f(m,k,:) = U(m,k) * H(:,k,m,ite) * G(m,k)';
                    v(m,k) = log2( U(m,k) );
                end
            end

        %%%%%%%%%%% CVX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cvx_begin quiet
            variable P_opt(Nt, K) complex;

            expression S;
            S = 0;
            for k = 1:K
                p_Psi = 0;
                Psi_bar = reshape(mean(Psi(:,k,:,:),1),Nt,Nt);
                Psi_bar = 0.5 * (Psi_bar + Psi_bar'); %to make it H
                for i = 1:K
                    p_Psi = p_Psi + ( P_opt(:,i)' * Psi_bar * P_opt(:, i) );
                end
                temp_f = mean(f(:,k,:),1);
                S = S + p_Psi + ( sigma * mean( t(:, k) ) ) - ( 2 * real(reshape...
                    (temp_f,Nt,1)'*P_opt(:, k))) + mean( U(:, k) )...
                    - mean( v(:, k) );
            end

            expression pwr;
            pwr = 0;
            for k = 1:K
                pwr = pwr + sum_square_abs( P_opt(:, k) );
            end

            minimize S;
            subject to;
            pwr <= Pt;

            cvx_end
        %%%%%%%%%%% CVX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            SVal = [SVal, S]; %value of cost functions over all iterations
            n = n + 1;
            if ~exist('Sold','var')
                a = 1;
            elseif n > 5 %|| abs(Sold - S) < .01 %to control while loop
                break;
            end

            Sold = S;
            P = P_opt;
            %a = 0;
        end
        %plot(1:n, SVal,'-b*');
        gamma = zeros(K,1);
        for m = 1:MV
            
            for k = 1:K
                Ik = 0;
                for i = 1:K
                    if i ~= k
                        Ik = Ik + abs( HVal(:,k,m)' * P(:,i) )^2 + sigma; %Interfearence power
                    end
                end
                gamma(k) = abs( HVal(:,k,m)' * P(:,k) )^2 / Ik;
                Rate(m,k) = log2(1+gamma(k));
            end
        end
        SRate(ite, :) = mean(Rate,1); %sum rate for one set of error samples
    end
    AvRate(p) = sum(mean(SRate,1)); %averaged sum rate
end
plot(SNR, AvRate, '-b*');