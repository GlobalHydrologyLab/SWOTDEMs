%% profConstrain.m 
% setup and solve slope-constrained lsqlin equations.

clear
close all


load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacDataV3.mat')
zField = 'geoHeight';
% simulated = trimFields(simulated,15:30); %test small area

simAllign = nodeAllign(simulated);
truthAllign = nodeAllign(truth);

nX = size(simAllign.(zField),1);
nDays= size(simAllign.(zField),2);
% nObs = size(simAllign.(zField),2);
nObs = sum(sum(~isnan(simAllign.(zField))));

%% do eet

%temporary placeholder for std
%in this case, more obs is better so no inverse. less is better with std,
%so make sure to add 1/std when replacing.
simAllign.nObsScaled = (simAllign.nObs - min(min(simAllign.nObs))) / max(max(simAllign.nObs));


C = zeros(nObs,nX); %weights
d = zeros(nObs,1);
A = zeros(nObs, nX);
b = zeros(nObs,1);

lastRow = 1;
for i = 1:nX
    ix = ~isnan(simAllign.nObsScaled(i,:));
    nTmp = sum(ix);
    endRow = lastRow+nTmp-1;
%     C(lastRow:endRow,i) = simAllign.nObsScaled(i,ix); %weight matrix.
%     d(lastRow:endRow,1) = simAllign.(zField)(i,ix) .* simAllign.nObsScaled(i,ix);
    C(lastRow:endRow,i) = zeros(endRow-lastRow + 1,1) + 1;  
    d(lastRow:endRow,1) = simAllign.(zField)(i,ix);
    
    if i < nX
        A(lastRow:endRow,i) = -1;
        A(lastRow:endRow,i+1) = 1;
    end
    lastRow = endRow + 1;
end

%%

tic
z = lsqlin(C,d,A,b);
toc

meanZ = nanmean(simAllign.(zField),2); %simple average
truthZ = nanmean(truthAllign.(zField),2);

%% plots

figure
plot(meanZ);
hold on
plot(z,'r','Linewidth',2);
plot(truthZ,'k','Linewidth',2)
hold off