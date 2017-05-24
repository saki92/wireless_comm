clear;

%% Declaration
P = 4; %No of transmit antennas
K = 4; %No of users (receivers)
Q = 2; %No of antennas per user

Nch = 20; %No of channel realizations

%% SNR vs sum-rate plot

%Noise covariance%
%R = eye(Q);
Etx = 1; %transmit power


SNR = (-10:5:30);
Etxm = 10.^(SNR/10);

for s = 1:length(SNR)
    Etx = Etxm(s);
    for ch = 1:Nch
        %Channel Matrix%
        Hm = ( 1/sqrt(2) ) * ( randn(Q*K,P) + 1i*randn(Q*K,P) );
        %Init transmit matrix%
        B = 1/sumsqr(abs(Hm'))*Hm'; %zero forcing matched filter
        H = mat2cell(Hm,Q*ones(1,K)); %convert matrix to cell for ease of iteration
        for n = 1:30 %fixed no of iterations
            [B,Rate] = sumRateCompute(B,H,Hm,P,K,Q,Etx);
            %calc no of iterations
            n = n + 1;
            sumRate(n) = abs(sum(Rate));
        end
        s_sumRate(ch) = sumRate(30);
    end
    f_sumRate(s) = mean(s_sumRate);
end
plot(SNR,f_sumRate,'-rx');
xlabel('SNR (dB)');
ylabel('Sum Rate');
grid on;