% SVD_multiSect_varReach.m
% Extending idea of SVD application to reduce errors in each profile
% individually to work for full river profiles instead of subsections.
% Reaches are defined by equiReach.m, which takes in total number of nodes
% and desired number of nodes per reach.

%--------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------
% X look into removing erroneous data from orbit 527 before SVD. The large
%       (and obvious) errors can dominate the decomposition results and
%       make comparisons of errors kind of disingenuous.
% X-- Hard-coded removal for now. I think this is fine. 
%
% - make sure # of singular values chosen is providing best results.
%
% - General cleaning of code.
%
% - reimagine error summary plots. Bar graph is good up to about ~15 groups
%        on the x-axis, then it becomes difficult to read.
%--------------------------------------------------------------------------

clear
close all

% % Sac
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacDataV3.mat')
zField = 'geoHeight';
% hard-coded removal of far range data from pass 527
for i = 10:17
    simulated(i).geoHeight(339:487) = NaN;
end
simulated(6) = [];
truth(6) = [];

% % Po
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoData.mat')
% zField = 'nHeight';
% % % % hard-coded removal of far range data
% simulated = trimFields(simulated,1:400);
% truth = trimFields(truth,1:400);

% % PoV2
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoDataV2.mat')
% zField = 'nHeight';
% % % % hard-coded removal of far range data
% simulated(18:35) = trimFields(simulated(18:35),1:400);
% truth(18:35) = trimFields(truth(18:35),1:400);

% % Tanana
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/transformedTananaData.mat')
% zField = 'nHeight';

clearvars -except simulated truth zField

%Group data into matrices without gaps, intelligently deleting data so that
%all sections have >= sectMin rows
nProf = length(simulated);
sectMin = 25;
simAllign = nodeAllign(simulated);

nodeRng = [min(simAllign.node(1,:)), max(simAllign.node(end,:))];
truth = trimFields(truth,nodeRng);

%--------------------------------------------------------------------------
% Create reaches
%--------------------------------------------------------------------------
nodePerReach = 25;
reachVec = equiReach(range(nodeRng)+1,nodePerReach);
[simulated.reach] = deal(reachVec);
[truth.reach] = deal(reachVec);
%--------------------------------------------------------------------------

truthAllign = nodeAllign(truth);
simAllign.(zField)(simAllign.(zField)==-9999) = NaN; %missing data value
truthAllign.(zField)(truthAllign.(zField)==-9999) = NaN; %missing data value

[section, zArray] = subsectByObs(simAllign.(zField),sectMin);
observedBy = ~isnan(zArray);


%init. matrices for storing section data
dim = size(simAllign.sCoord);
z2All = nan(dim);
sAll = nan(dim);

for r = min(section):max(section)
    inSect = section == r;
    
    simReach = trimFields(simAllign,inSect);
    truthReach = trimFields(truthAllign,inSect);
    
    %assemble full rectangular matrices of s,z data 
    z = zArray(inSect,:);
    [z, delCol] = nanRows(z,2);
    s = simAllign.sCoord(inSect,~delCol);

    %replace any missing s-values with mean of that node.
    [mMiss, nMiss] = find(isnan(s));
    for i = 1:length(mMiss)
        s(mMiss(i),nMiss(i)) = nanmean(s(mMiss(i),:));
    end
    
    % fit polynomial to all data in reach
    nProfSect = size(z,2);
    polyOrder = 1;
    sv = reshape(s,[],1);
    zv = reshape(z,[],1);
    [m,~,mu] = polyfit(sv,zv,polyOrder);
    
    %SVD on residuals
    zhat = polyval(m,s,[],mu);
    zresid = reshape(z-zhat, [], nProfSect);
    
    [U,S,V] = svd(zresid);
    
    %----------------------------------------------------------------------
    % Determine best rank 
    %
    % This section needs to be finished- both break in magnitude and IQR
    % approach perform better than the other in some cases. Maybe some kind
    % of consensus with a third method could work.
    %----------------------------------------------------------------------
    SVal = diag(S);
    
    % WORKING IDEAS:
    
    %largest decrease in magnitude
%     [~,iSV] = min(diff(S));

    %inter quartile range approach.
    IQR = iqr(SVal);
    try
        iSV(r) = find(SVal >= median(SVal) + 2*IQR,1,'last');
    catch
        iSV(r) = 1;
    end
    
    %cumulative percent of sing. val. threshold.
%     p = 0.5;
%     iSV = find(cumsum(S)./sum(S) >= p,1,'first')
    %----------------------------------------------------------------------

    %now modify S, removing smaller components.
    S2 = S;
    S2(iSV+1:end,:) = 0;
    
    %recombine
    z2resid = U*S2*V';
    z2 = z2resid + polyval(m,s,[],mu);
  
    % join section data for later comparison
    z2All(inSect,~delCol) = z2;
    sAll(inSect,~delCol) = s;
end

skm = nanmean(sAll,2)/1000;
zAll = simAllign.(zField);

%--------------------------------------------------------------------------
% Reach Stats
%--------------------------------------------------------------------------
reaches = unique(reachVec)';

simStats.slopeErr = nan(length(reaches),nProf);
svdStats.slopeErr = nan(length(reaches),nProf);
simStats.reachAvgZ = nan(length(reaches),nProf);
svdStats.reachAvgZ = nan(length(reaches),nProf);
reachAvgTruth = nan(length(reaches),nProf);

for r = reaches
    RL(r) = range(skm(reachVec==r));
    for p = 1:nProf
        inReach = truthAllign.reach(:,p) == r & ~isnan(z2All(:,p));
        if sum(inReach)>=2
            fitSim = polyfit(sAll(inReach,p),zAll(inReach,p),1);
            fitSVD = polyfit(sAll(inReach,p),z2All(inReach,p),1);

            inReach = truth(p).reach == r;
            fitTruth = polyfit(truthAllign.sCoord(inReach,p), ... 
                truthAllign.(zField)(inReach,p),1);

            simStats.slopeErr(r,p) = fitTruth(1) - fitSim(1);
            svdStats.slopeErr(r,p) = fitTruth(1) - fitSVD(1);
            simStats.relSlopeErr(r,p) = simStats.slopeErr(r,p) / fitTruth(1);
            svdStats.relSlopeErr(r,p) = svdStats.slopeErr(r,p) / fitTruth(1);
            
            simStats.reachAvgZ(r,p) = nanmean(zAll(inReach,p));
            svdStats.reachAvgZ(r,p) = nanmean(z2All(inReach,p));
            reachAvgTruth(r,p) = nanmean(truthAllign.(zField)(inReach,p));
        end
    end
end

simStats.reachZErr = reachAvgTruth - simStats.reachAvgZ;
svdStats.reachZErr = reachAvgTruth - svdStats.reachAvgZ;

%% 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

%rothko section plot
% imagesc(observedBy .* section)
% xlabel('Profile')
% ylabel('Node')
% c = lines;
% c = c(1:max(section),:);
% colormap([1 1 1; c])

%singular values
figure()
bar(diag(S))
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
title('Singular Values of Elevation Residuals')

%original and approx profiles
handle = figure();
subplot(2,1,1);
plot(skm,zArray)
hold on
plot(skm,truth(3).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Original Simulated Profiles')

subplot(2,1,2);
plot(skm,z2All)
hold on
plot(skm,truth(3).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Low-Rank Approximation Profiles')

for r = reaches
    ir = find(reachVec == r,1,'last');
    xr = skm(ir);
    yr = truthAllign.(zField)(ir,3);
    
    subplot(2,1,1)
    plot([xr xr],[yr yr], 'r*','LineWidth',2)
    
    subplot(2,1,2)
    plot([xr xr],[yr yr], 'r*','LineWidth',2)
end

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position', [0.2, 0.1, 0.6, 0.9])


zErr = zAll - truthAllign.(zField);
z2Err = z2All - truthAllign.(zField);

%epdf of all node errors
figure()
ksdensity(reshape(zErr,[],1))
hold on
ksdensity(reshape(z2Err,[],1))
xlabel('Elevation Error (m)')
title('Empirical PDF of Node Errors')
legend('Original Data','Low Rank')

% reach slope errors
simStats.slopeMAE = nanmean(abs(simStats.slopeErr),1);
svdStats.slopeMAE = nanmean(abs(svdStats.slopeErr),1);

simStats.slopeRMSE = sqrt(nanmean(simStats.slopeErr.^2,1));
svdStats.slopeRMSE = sqrt(nanmean(svdStats.slopeErr.^2,1));

simStats.slopeRRMSE = sqrt(nanmean(simStats.relSlopeErr.^2,1)).*100;
svdStats.slopeRRMSE = sqrt(nanmean(svdStats.relSlopeErr.^2,1)).*100;

figure()
plot(simStats.slopeRRMSE', svdStats.slopeRRMSE','k.','MarkerSize',20)
% plot(abs(simStats.relSlopeErr), abs(svdStats.relSlopeErr),'ko')
grid on
axis square
axLim = get(gca,'XLim');
axLim(1) = 0;
hold on
plot(axLim,axLim,'k-')
set(gca,'XLim', axLim); set(gca,'YLim',axLim);
title('Reach-level slope errors')
xlabel('Simulated RRMSE(%)')
ylabel('LRA RRMSE(%)')

% reach elev. errors
simStats.reachZMAE = nanmean(abs(simStats.reachZErr),1);
svdStats.reachZMAE = nanmean(abs(svdStats.reachZErr),1);

simStats.reachZRMSE = sqrt(nanmean(simStats.reachZErr.^2,1));
svdStats.reachZRMSE = sqrt(nanmean(svdStats.reachZErr.^2,1));

figure()
plot(simStats.reachZRMSE', svdStats.reachZRMSE','k.','MarkerSize',20)
% plot(abs(simStats.reachZErr), abs(svdStats.reachZErr),'ko')
grid on
axis square
axLim = get(gca,'YLim');
% axLim(1) = 0;
hold on
plot(axLim,axLim,'k-')
set(gca,'XLim', axLim); set(gca,'YLim',axLim);
title('Reach-level elevation errors')
xlabel('Simulated (m)')
ylabel('LRA (m)')

%--------------------------------------------------------------------------


