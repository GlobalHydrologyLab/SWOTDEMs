% testSVD.m
% trying singular value decomposition method.

%   TO DO:
% - extend to loop through sections as defined by subsectByObs.m
% x- join those sections (added to subsectByObs.m)
% x- might need to group shorter sections and remove extra observations
% - deals with range effects
% - figure out best number of singular values to use.
% - evaluate smoothing methods/windows

clear
close all

% % Sac
% load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Sacramento/transformedSacDataV2.mat')
% zField = 'geoHeight';
% for i = 1:length(truth)
%     %true where widths are missing or height is default -9999 value.
%     iBadHeight = isnan(truth(i).nWidth) | truth(i).geoHeight == -9999; 
%     
%     simulated(i).reach(iBadHeight) = NaN;
%     simulated(i).node(iBadHeight) = NaN;
%     simulated(i).easting(iBadHeight) = NaN;
%     simulated(i).northing(iBadHeight) = NaN;
%     simulated(i).nHeight(iBadHeight) = NaN;
%     simulated(i).nWidth(iBadHeight) = NaN;
%     simulated(i).geoHeight(iBadHeight) = NaN;
%     simulated(i).sCoord(iBadHeight) = NaN; 
%     simulated(i).nCoord(iBadHeight) = NaN;
%     
%     
%     truth(i).reach(iBadHeight) = NaN;
%     truth(i).node(iBadHeight) = NaN;
%     truth(i).easting(iBadHeight) = NaN;
%     truth(i).northing(iBadHeight) = NaN;
%     truth(i).nHeight(iBadHeight) = NaN;
%     truth(i).nWidth(iBadHeight) = NaN;
%     truth(i).geoHeight(iBadHeight) = NaN;
%     truth(i).sCoord(iBadHeight) = NaN; 
%     truth(i).nCoord(iBadHeight) = NaN;
% end

% % Tanana
% load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Tanana/transformedTananaData.mat')
% zField = 'nHeight';

% % % Po
load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Po/transformedPoDataV2.mat')
zField = 'nHeight';

clearvars -except simulated truth zField
alligned = nodeAllign(simulated);

%% Sectify!
sectMin = 25;
[section, zArray] = subsectByObs(alligned.(zField),sectMin);
observedBy = ~isnan(zArray);

imagesc(observedBy .* section)
xlabel('Profile')
ylabel('Node')
c = lines;
c = c(1:max(section),:);
colormap([1 1 1; c])

%% SVD
z2Avg = []; %concatenation of processed sections.
zChk = [];

handle = figure();
subplot(2,1,1); hold on
subplot(2,1,2); hold on

for i = min(section):max(section)
    allignedIDs = section == i;
    sectNodes = alligned.nodeVec(allignedIDs);
    
    z = zArray(allignedIDs,:);
    z = nanRows(z,2); %delete 'NaN space' to form nice matrices 
    
    
    
    [U,S,V] = svd(z,0);

    %zChk should be the same as z.
    zChkSect = U*S*V';
    zChk = [zChk; mean(zChkSect,2)];

    %truncate S, removing smaller components.
    k=2;
    S2 = S;
    S2(k+1:end,:) = 0;
    
    
%     R=S2*V';
%     R=R(1:k,:);
%     q = iqr(R);
%     m = median(R);
%     iqrOutlier = R < (m - 1.5*q) | R > (m + 1.5*q);
%     V(iqrOutlier,:) = NaN;

    
    zSect = U*S2*V';
    
    zSectAvg = nanmedian(zSect,2); 
    z2Avg = [z2Avg; zSectAvg];
    
    subplot(2,1,1)
    plot(sectNodes, z, 'color', c(i,:));

end

%trim sim and truth, take find avg.
svdNodes = ~isnan(section) .* alligned.nodeVec;
svdNodes(svdNodes==0) = [];

truth = trimFields(truth, svdNodes);
truthAvg = nodeAvg3_1(truth, zField);
simulated = trimFields(simulated, svdNodes);
simAvg = nodeAvg3_1(simulated, zField);

%avg bias
biasAvg = nanmean(simAvg.(zField) - truthAvg.(zField));
simAvg.(zField) = simAvg.(zField) - biasAvg;

%smooth svd value
svdSmooth = smooth(z2Avg,21,'rloess');

%svdSmooth bias
biasSvdSmooth = nanmean(svdSmooth - truthAvg.(zField));
svdSmooth = svdSmooth - biasSvdSmooth;

%svd bias
biasSvd = nanmean(z2Avg - truthAvg.(zField));
z2Avg = z2Avg - biasSvd;

subplot(2,1,1)
plot(svdNodes,z2Avg,'k','Linewidth',2')
plot(svdNodes,truthAvg.(zField), 'b', 'Linewidth',2)

subplot(2,1,2)
plot(svdNodes,z2Avg,'k','Linewidth',2')
plot(svdNodes,truthAvg.(zField), 'b', 'Linewidth',2)
plot(svdNodes,svdSmooth,'r','Linewidth',2)
legend('Truncated SVD', 'Averaged Input', 'Smoothed TSVD','Location','NorthWest')
% skm = simAvg.sCoord/1000;
% 
% handle = figure();
% subplot(2,1,1);
% plot(skm,zChk)
% hold on
% % plot(skm,truthAvg.(zField),'k','Linewidth',2)
% 
% subplot(2,1,2);
% plot(skm,z2Avg)
% hold on
% % plot(skm,truthAvg.(zField),'k','Linewidth',2)
% % plot(simAvg.(zField),'ro')
% plot(skm,svdAvg,'r-','Linewidth',2)

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position', [0.2, 0.1, 0.6, 0.9])

% errors
RMSEavg = sqrt(mean((simAvg.(zField) - truthAvg.(zField)).^2));
MAEavg = mean(abs((simAvg.(zField) - truthAvg.(zField)).^2));

RMSEsvd = sqrt(mean((svdSmooth - truthAvg.(zField)).^2));
MAEsvd = mean(abs(svdSmooth - truthAvg.(zField)));
