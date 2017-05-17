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

n = 0;
a = 1;
while a == 1
    [B,Rate] = sumRateCompute(B,H,Hm,P,K,Q,Etx);
    n = n + 1;
    if ~exist('Bold','var')
        a = 1;
    elseif n > 30 %|| abs(sumsqr(abs(Bold))-sumsqr(abs(B))) < .000001 %to control while loop
        a = 0;
    end

    Bold = B;
    sumRate(n) = abs(sum(Rate));
end
plot([0:n-1],sumRate,'-kx');

%% SNR vs sum-rate plot
SNR = (-10:5:30);
Etxm = 10.^(SNR/10);
H = mat2cell(Hm,Q*ones(1,K)); %convert matrix to cell for ease of iteration
for s = 1:length(SNR)
    Etx = Etxm(s);
    for n = 1:30 %fixed no of iterations
        [B,Rate] = sumRateCompute(B,H,Hm,P,K,Q,Etx);
        %calc no of iterations
        n = n + 1;
        sumRate(n) = abs(sum(Rate));
    end
    f_sumRate(s) = sumRate(30);
end
plot(SNR,f_sumRate,'-rx');