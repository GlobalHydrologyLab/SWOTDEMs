function [d] = SVRecomp(U,S,V,k)
%SVRecomp.m
%takes in U,S,V components of singular value decomposition and reassembles
%them using only factors in vector k. 

Sk = zeros(size(S));
Sk(k,k) = S(k,k);

d = U * Sk * V';
end

