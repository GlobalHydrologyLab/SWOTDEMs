%% profConstrain.m 
% setup and solve slope-constrained lsqlin equations.

clear
close all


load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacDataV4.mat')
zField = 'geoHeight';

% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPo_3Pass.mat')
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/transformedTananaData.mat')
% zField = 'nHeight';
% simulated = trimFields(simulated,15:30); %test small area

simAllign = nodeAllign(simulated);
truthAllign = nodeAllign(truth);
meanZ = nanmean(simAllign.(zField),2); %simple average
truthZ = nanmean(truthAllign.(zField),2);

nX = size(simAllign.(zField),1);
% nDays= size(simAllign.(zField),2);
nObs = sum(sum(~isnan(simAllign.(zField))));

%% set up

%temporary placeholder for std
%in this case, more obs is better so no inverse. less is better with std,
%so make sure to do 1/std when replacing.
simAllign.nObsScaled = (simAllign.nObs - min(min(simAllign.nObs))) / max(max(simAllign.nObs));


G = zeros(nObs,nX); %weights
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
    G(lastRow:endRow,i) = zeros(endRow-lastRow + 1,1) + 1;  
    d(lastRow:endRow,1) = simAllign.(zField)(i,ix);
    
    if i <= nX - 1
        A(lastRow:endRow,i) = -1;
        A(lastRow:endRow,i+1) = 1;
    end
%     if i <= nX - 2
%        A(lastRow:endRow,i+1) = A(lastRow:endRow,i+1) + 1;
%        A(lastRow:endRow,i+2) = -1; 
%     end
    lastRow = endRow + 1;
end

%% do it
% tic
% z = lsqlin(C,d,A,b);
% toc

[U,S,V] = svd(G,'econ');
Sinv = diag(1./diag(S));
tic
zChk = V * Sinv * U' * d;
toc

% % note:
% % G = U*C*X'
% % A = V*S*X';


% Should A constrain data or model?
A = [];
for i = 1:nX-1
A(i,i) = -1;
A(i,i+1) = 1;
end
A(end+1,:) = 0;

%%
tic
[U,V,X,C,S] = gsvd(G,A);
toc



%% plots


figure
plot(meanZ);
hold on
plot(z,'r','Linewidth',2);
plot(truthZ,'k','Linewidth',2)
hold off