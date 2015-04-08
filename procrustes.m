function W = procrustes(A, B)
% Find W such that A*W approximates B.
    [U, S, V] = svd(A'*B);
    W = U*V';
end