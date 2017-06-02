clear;
clc;

Nt = 2; %no of transmit antennas
K = 2; %no of receivers
Q = 1; %no of antennas per receiver
M = 10; %no of channel realizations
sigma = 0.5; %noise covariance

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
        eta{m+1,k} = t(m+1,k) * Hm(m*Nt+1:(m*Nt+1)+Nt-1,k) * ...
            Hm(m*Nt+1:(m*Nt+1)+Nt-1,k)';
    end
end
