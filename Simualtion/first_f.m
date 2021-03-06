clear;

%% Declaration
P = 4; %No of transmit antennas
K = 2; %No of users (receivers)
Q = 2; %No of antennas per user

%Channel Matrix%
% Hm = zeros(Q*K,P);
% for k = 1:K*Q
%     Hm(k,:) = randn(1,P) + 1i*randn(1,P); %vertical divided
% end

Hm = randn(Q*K,P) + 1i*randn(Q*K,P);

%Init transmit matrix%
%B = ones(P,K*Q); %all ones
B = 1/sumsqr(abs(Hm'))*Hm'; %zero forcing matched filter

%Noise covariance%
%R = eye(Q);
Etx = 1; %transmit power

%% Compute WMMSE matrix

H = mat2cell(Hm,Q*ones(1,K)); %convert matrix to cell for ease of iteration

n = 1;
a = 1;
while a == 1
    [B,Rate] = sumRateCompute(B,H,Hm,P,K,Q,Etx);
    sumRate(n) = abs(sum(Rate));
    if n == 1
        a = 1;
    elseif n > 300 || abs(sumsqr(abs(sumRate(n)))-sumsqr(sumRate(n-1))) < 10^-22 %to control while loop
        a = 0;
    end

    n = n + 1;
    
end
plot([1:length(sumRate)],sumRate,'-kx');