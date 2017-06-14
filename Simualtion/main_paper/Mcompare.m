M = [1,5,10,15,20,30];
SNR = 5:3:35;
hold on,
for m = 1:length(M)
    RateVec{m} = sumRate_2_2(m);
    plot(SNR, AvRate, '-*');
    m
end
hold off;
