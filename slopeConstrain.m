function [x] = slopeConstrain(d,maxSlope)
% x = slopeConstrain(d,maxSlope)
%
%Takes in data vector d, sets up and solves lsqlin() such that:
%
%d(i) - d(i+1) <= maxSlope
%
%maxSlope has a default value of 0.
%
%Read documentation for lsqlin for more info.

if nargin < 2 
    maxSlope = 0;
end

d=d(:); %force column vector.
n = numel(d);

x = nan(n,1);

[d, delId] = nanRows(d);
n = n-sum(delId);

C = diag(ones(n,1));
b = zeros(n,1) + maxSlope;
A = zeros(n,n);

for j = 1:n-1
    A(j,j) = -1;
    A(j,j+1) = 1;
end

options = optimset('display','off');
[xCol,~,~,exitFlag] = lsqlin(C,d,A,b,[],[],[],[],[],options);

if exitFlag ~= 1
    warning(['lsqlin exit flag was ' num2str(exitFlag) ...
        '. Check lsqlin documentation for cause of this flag.'])
end

x(~delId) = xCol;

end

