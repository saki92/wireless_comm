clear;
clc;
Nt = 2; %no of transmit antennas
K = 2; %no of receivers
M = [1,10,100,1000];
SNR = 5:3:35;
noofit = 1;
alpha = 0.6;
for av = 1:100
    Hcap = ( 1/sqrt(2) ) * ( randn(Nt,K) + 1i*randn(Nt,K) );
    % P = ( 1 / sumsqr(abs(reshape(Hcap,2,2)')) * reshape(Hcap,2,2)' ).'; %initial Tx vectors for all samples
    % P = ones(Nt,K);
    P = conj(Hcap);

    for m = 1:length(M)
        RateVec(:,m,av) = sumRate_2_2(M(m),SNR,Hcap,noofit,P,alpha);
    end
    
%     RateVec_M1(:,1,av) = sumRateM1(1,SNR,Hcap,noofit,P,alpha,1);
%     RateVec_M1(:,2,av) = sumRateM1(1,SNR,Hcap,noofit,P,alpha,0);

    av
end
RateVec = mean(RateVec,3);
%% plot
plotStyle = {'-*b','-xk','-or','-dg'};
legendinfo = {'M = 1','M = 10','M = 100','M = 1000'};
hold on,
for m = 1:length(M)
    plot(SNR, RateVec(:,m), plotStyle{m}');
end
% plot(SNR, mean(RateVec_M1(:,2),3), ':*k');
% plot(SNR, mean(RateVec_M1(:,1),3), ':*c');
grid on;
legend(legendinfo);
hold off;