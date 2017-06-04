clear;
clc;

Nt = 2; %no of transmit antennas
K = 2; %no of receivers
Q = 1; %no of antennas per receiver
M = 10; %no of channel realizations
sigma = 0.5; %noise covariance
Pt = 0.5; %total transmit power

Hm = [];
for m = 1:M
    Hm = [ Hm ; ( 1/sqrt(2) ) * ( randn(Nt,Q*K) + 1i*randn(Nt,Q*K) ) ];
end

P = ( 1 / sumsqr(abs(Hm')) * Hm' ).';

T = [];
for k = 1:K
    for m = 0:M-1
        for i = 1:K
            if i ~= k
                Ik = abs( Hm(m*Nt+1:(m*Nt+1)+Nt-1,k)' * ...
                    P(m*Nt+1:(m*Nt+1)+Nt-1,i) )^2 + sigma; %Interfearence power
            end
        end
        T(m+1,k) = abs( Hm(m*Nt+1:(m*Nt+1)+Nt-1,k)' * ...
            P(m*Nt+1:(m*Nt+1)+Nt-1,k) )^2; %Average power
        G(m+1,k) = P(m*Nt+1:(m*Nt+1)+Nt-1,k)' *...
            Hm(m*Nt+1:(m*Nt+1)+Nt-1,k) / T(m+1,k); %Equalizer matrix
        U(m+1,k) = inv( T(m+1,k) \ Ik); %MMSE weights
        t(m+1,k) = U(m+1,k) * abs( G(m+1,k) )^2;
        eta(m+1,k,:,:) = t(m+1,k) * Hm(m*Nt+1:(m*Nt+1)+Nt-1,k) * ...
            Hm(m*Nt+1:(m*Nt+1)+Nt-1,k)';
        f(m+1,k,:) = U(m+1,k) * Hm(m*Nt+1:(m*Nt+1)+Nt-1,k) * G(m+1,k)';
        v(m+1,k) = log2( U(m+1,k) );
    end
end

cvx_begin quiet
variable P_opt(Nt, K) complex;

expression S;
for k = 1:K
    p_eta = 0;
    for i = 1:K
        p_eta = p_eta + P_opt(:,i)' * mean( eta(:, k) ) * P_opt(:, i);
    end
    temp_f = mean(f(:,k,:),1);
    S = S + p_eta + ( sigma * mean( t(:, k) ) ) - ( 2 * real(reshape(temp_f(1,1,:),2,1)'...
        *P_opt(:, k))) + mean( U(:, k) ) - mean( v(:, k) );
end

expression pwr;
for k = 1:K
    pwr = pwr + sum_square_abs( P_opt(:, k) );
end

minimize S;
subject to
pwr <= Pt;

cvx_end
