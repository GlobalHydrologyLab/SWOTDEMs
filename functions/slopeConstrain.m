function [x] = slopeConstrain(dIn,maxDiff,weight)
% x = slopeConstrain(d,maxDiff)
%
%Takes in data d, sets up and solves lsqlin() such that:
%
%x(i) - x(i+1) <= maxDiff
%
%d can be a matrix or a column vector. rows are considered to be the same x
%value in the regression. Calling slopeConstrain on a matrix will use all
%non-NaN values in the row to estimate the row value. maxDiff has a
%default value of 0. 
%
%Read documentation for lsqlin for more info.


if ~exist('maxDiff','var')
    maxDiff = 0;
end

if ~exist('weight','var')
    weight = ones(size(dIn));
else
    %set any 0 or NaN weights to 1.
    weight(weight==0) = 1;
    weight(isnan(weight)) = 1;
end

observed = ~isnan(dIn);
missingNodes = ~any(observed,2);
nObs = sum(sum(observed));

nOrigNodes = size(dIn,1);
dIn(missingNodes,:) = [];
nNodes = size(dIn,1);
observed = ~isnan(dIn);

C = zeros(nObs,nNodes);
d = zeros(nObs,1);
A = zeros(nNodes);
b = zeros(nNodes,1) + maxDiff;


% Objective function:   ||Cx - d||2
% contrained by:        Ax <= b

firstRow = 1;
for i = 1:nNodes
    nd = sum(~isnan(dIn(i,:)));
    
    if nd ~=0
        endRow = firstRow+nd-1;

        C(firstRow:endRow,i) = ones(nd,1) ./ weight(i,observed(i,:))';
        d(firstRow:endRow,1) = dIn(i,observed(i,:)) ./ weight(i,observed(i,:));
        
        firstRow = endRow + 1;

        %first derivative matrix
        if i < nNodes
            A(i,i+1) = 1;
            A(i,i) = -1;
        end  
        
    end
end

options = optimset('display','off');
[xSolve,~,~,exitFlag] = lsqlin(C,d,A,b,[],[],[],[],[],options);

if exitFlag ~= 1
    warning(['lsqlin exit flag was ' num2str(exitFlag) ...
        '. Check lsqlin documentation for cause of this flag.'])
end

x = nan(nOrigNodes,1);
x(~missingNodes) = xSolve;

end

