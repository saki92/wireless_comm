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
    B = mat2cell(B,P,Q*ones(1,K)); %convert matrix to cell for ease of iteration
    A = cell(1,K);
    W = cell(1,K);
    
    %finding noise covariance matrix%
    for k = 1:K
        t = zeros(Q);
        for i = 1:K
            if i~=k
                t = t + (H{k}*B{i}*B{i}'*H{k}');
            end
        end
        R{k} = eye(Q) + t;
    end
    
    for k = 1:K
        %MMSE receive matrix
        A{k} = B{k}'*H{k}' / ((H{k}*B{k}*B{k}'*H{k}' + R{k}));
        %MSE weights
        E = inv(eye(Q) + B{k}'*H{k}'/R{k}*H{k}*B{k});
        W{k} = eye(Q)/E;
        Rate(k) = log(det(inv(E)));
    end
    
    A = cell2blk(A); %convert cell to block diag
    W = cell2blk(W); %convert cell to block diag
    %MMSE transmit matrix
    B = (Hm'*(A')*W*A*Hm + (trace(W*A*(A'))/Etx)*eye(P))\Hm'*A'*W;
    b = sqrt(Etx/trace(B*B'));
    B = b*B;
    %calc no of iterations
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
        B = mat2cell(B,P,Q*ones(1,K)); %convert matrix to cell for ease of iteration
        A = cell(1,K);
        W = cell(1,K);

        %finding noise covariance matrix%
        for k = 1:K
            t = zeros(Q);
            for i = 1:K
                if i~=k
                    t = t + (H{k}*B{i}*B{i}'*H{k}');
                end
            end
            R{k} = eye(Q) + t;
        end

        for k = 1:K
            %MMSE receive matrix
            A{k} = B{k}'*H{k}' / ((H{k}*B{k}*B{k}'*H{k}' + R{k}));
            %MSE weights
            E = inv(eye(Q) + B{k}'*H{k}'/R{k}*H{k}*B{k});
            W{k} = eye(Q)/E;
            Rate(k) = log(det(inv(E)));
        end

        A = cell2blk(A); %convert cell to block diag
        W = cell2blk(W); %convert cell to block diag
        %MMSE transmit matrix
        B = (Hm'*(A')*W*A*Hm + (trace(W*A*(A'))/Etx)*eye(P))\Hm'*A'*W;
        b = sqrt(Etx/trace(B*B'));
        B = b*B;
        %calc no of iterations
        n = n + 1;
        if ~exist('Bold','var')
            a = 1;
        elseif n > 30 %|| abs(sumsqr(abs(Bold))-sumsqr(abs(B))) < .000001 %to control while loop
            a = 0;
        end

        Bold = B;
        sumRate(n) = abs(sum(Rate));
    end
    f_sumRate(s) = sumRate(30);
end
plot(SNR,f_sumRate,'-rx');