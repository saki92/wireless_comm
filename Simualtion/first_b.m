clear;

%% Declaration
P = 4; %No of transmit antennas
K = 2; %No of users (receivers)
Q = 2; %No of antennas per user

%Channel Matrix%
Hm = zeros(Q*K,P);
for k = 1:K*Q
    Hm(k,:) = randn(1,P); %vertical divided
end

%Init beamforming weights%
B = ones(P,K*Q); %horizontal divided

%Noise covariance%
R = eye(Q);
Etx = 1; %transmit power

%% Compute WMMSE matrix

H = mat2cell(Hm,Q*ones(1,K)); %convert matrix to cell for ease of iteration

n = 0;
a = 1;
while a == 1
    B = mat2cell(B,P,Q*ones(1,K)); %convert matrix to cell for ease of iteration
    A = cell(1,K);
    W = cell(1,K);
    for k = 1:K   
        %MMSE receive matrix
        A{k} = B{k}'*H{k}' / ((H{k}*B{k}*B{k}'*H{k}' + R));
        %MSE weights
        E = inv(eye(Q) + B{k}'*H{k}'/R*H{k}*B{k});
        W{k} = eye(Q)/E;
    end
    A = cell2blk(A); %convert cell to block diag
    W = cell2blk(W); %convert cell to block diag
    %MMSE transmit matrix
    B = (Hm'*(A')*W*A*Hm + (trace(W*A*(A'))/Etx)*eye(P))\Hm'*A'*W;
    %calc no of iterations
    n = n + 1;
    if ~exist('Bold','var')
        a = 1;
    elseif Bold ~= B & n < 400 %to control while loop
        a = 1;
    else
        a = 0;
    end
    Bold = B;
end