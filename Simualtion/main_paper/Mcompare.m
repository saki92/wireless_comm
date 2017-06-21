clear;
clc;
Nt = 2; %no of transmit antennas
K = 2; %no of receivers
M = [1,10,100,1000];
SNR = 5:3:35;
noofit = 1;
for av = 1:30
    Hcap = ( 1/sqrt(2) ) * ( randn(Nt,K) + 1i*randn(Nt,K) );
    % P = ( 1 / sumsqr(abs(reshape(Hcap,2,2)')) * reshape(Hcap,2,2)' ).'; %initial Tx vectors for all samples
    % P = ones(Nt,K);
    P = conj(Hcap);

    for m = 1:length(M)
        RateVec(:,m,av) = sumRate_2_2(M(m),SNR,Hcap,noofit,P);
    end
    av
end
RateVec = mean(RateVec,3);
%% plot
plotStyle = {'-*b','-xk','-or','-dg'};
legendinfo{'M = 1','M = 10','M = 100','M = 1000'};
hold on,
for m = 1:length(M)
    plot(SNR, RateVec(:,m), plotStyle{m}');
end
grid on;
legend(legendinfo);
hold off;