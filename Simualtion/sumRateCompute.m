function [ B,Rate ] = sumRateCompute( B,H,Hm,P,K,Q,Etx )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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

end

