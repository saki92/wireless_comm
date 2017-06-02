clear;
clc;

Nt = 2; %no of transmit antennas
K = 2; %no of receivers
Q = 1; %no of antennas per receiver
M = 10; %no of channel realizations

Hm = [];
for m = 1:M
    Hm = [ Hm ; ( 1/sqrt(2) ) * ( randn(Nt,Q*K) + 1i*randn(Nt,Q*K) ) ];
end

P = ( 1 / sumsqr(abs(Hm')) * Hm' ).';

T = [];
for k = 1:K
    for m = 0:M-1
        T(m+1,k) = abs( Hm(m*Nt+1:(m*Nt+1)+Nt-1,k)' * ...
            P(m*Nt+1:(m*Nt+1)+Nt-1,k) )^2; 
    end
end
