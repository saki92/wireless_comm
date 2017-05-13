function D = cell2blk(X)
if length(X) < 2
    error('Input cell size should be more than 2')
elseif length(X) == 2
    D = blkdiag(X{1},X{2});
else
    D = blkdiag(X{1},X{2});
    for n = 3:length(X)
        D = blkdiag(D,X{n});
    end
end