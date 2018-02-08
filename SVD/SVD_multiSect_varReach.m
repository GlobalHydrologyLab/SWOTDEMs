% SVD_multiSect_varReach.m
% Extending idea of SVD application to reduce errors in each profile
% individually to work for full river profiles instead of subsections.
% Reaches are defined by equiReach.m, which takes in total number of nodes
% and desired number of nodes per reach.

%--------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------
% - look into removing erroneous data from orbit 527 before SVD. The large
%       (and obvious) errors can dominate the decomposition results and
%       make comparisons of errors kind of disingenuous.
% X-- Hard-coded removal for now. I think this is fine. 
%
% - make sure # of singular values chosen is providing best results.
%
%--------------------------------------------------------------------------

clear
close all

targetRLkm = 10;

% river = 'Sac';
% river = 'Po';
% river = 'PoV2';
river = 'Tanana';

switch river
    case 'Sac'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacDataV4.mat')
        zField = 'geoHeight';
        % % hard-coded removal of far range data from pass 527
        for i = 10:17
            simulated(i).geoHeight(339:487) = NaN;
        end
        simulated(6) = [];
        truth(6) = [];

    case 'Po'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoData.mat')
        zField = 'nHeight';
        % % % hard-coded removal of far range data
        simulated = trimFields(simulated,1:400);
        truth = trimFields(truth,1:400);

    case 'PoV2'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoDataV2.mat')
        zField = 'nHeight';
        % % % hard-coded removal of far range data
        simulated(18:35) = trimFields(simulated(18:35),1:400);
        truth(18:35) = trimFields(truth(18:35),1:400);

    case 'Tanana'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/transformedTananaData.mat')
        zField = 'nHeight';
        simulated = trimFields(simulated,99:829);
        truth = trimFields(truth,99:829);
end


clearvars -except simulated truth zField targetRLkm


%Group data into matrices without gaps, intelligently deleting data so that
%all sections have >= sectMin rows
nProf = length(simulated);
sectMin = 25;
simAllign = nodeAllign(simulated);

nodeRng = [min(simAllign.node(1,:)), max(simAllign.node(end,:))];
truth = trimFields(truth,nodeRng);

truthAllign = nodeAllign(truth);
simAllign.(zField)(simAllign.(zField)==-9999) = NaN; %missing data value
truthAllign.(zField)(truthAllign.(zField)==-9999) = NaN; %missing data value

simAllign.sCoord = simAllign.sCoord;
truthAllign.sCoord = truthAllign.sCoord;

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

zErr = zAll - truthAllign.(zField);
z2Err = z2All - truthAllign.(zField);
%--------------------------------------------------------------------------
% Reach Stats
%--------------------------------------------------------------------------
%create reaches
reachVec = equiReach(skm,targetRLkm);
[simulated.reach] = deal(reachVec);
[truth.reach] = deal(reachVec);

reaches = unique(reachVec)';

%init stats matrices with NaNs.

simStats.slope = nan(length(reaches),nProf);
simStats.slopeErr = nan(length(reaches),nProf);
simStats.relSlopeErr = nan(length(reaches),nProf);
simStats.reachAvgZ = nan(length(reaches),nProf);

svdStats = simStats;
reachAvgTruth = nan(length(reaches),nProf);
truthSlope = nan(length(reaches),nProf);

for r = reaches
    RL(r) = range(skm(reachVec==r));
    for p = 1:nProf
        inReach = reachVec == r & ~isnan(z2All(:,p));
        if sum(inReach)>=2
            fitSim = polyfit(sAll(inReach,p),zAll(inReach,p),1);
            fitSVD = polyfit(sAll(inReach,p),z2All(inReach,p),1);

            inReach = truth(p).reach == r;
            fitTruth = polyfit(truthAllign.sCoord(inReach,p), ... 
                truthAllign.(zField)(inReach,p),1);
            
            simStats.slope(r,p) = fitSim(1) .*100000;
            svdStats.slope(r,p) = fitSVD(1) .*100000;
            truthSlope(r,p) = fitTruth(1) .*100000;

            simStats.slopeErr(r,p) = truthSlope(r,p) - simStats.slope(r,p);
            svdStats.slopeErr(r,p) = truthSlope(r,p) - svdStats.slope(r,p);
            simStats.relSlopeErr(r,p) = simStats.slopeErr(r,p) / truthSlope(r,p);
            svdStats.relSlopeErr(r,p) = svdStats.slopeErr(r,p) / truthSlope(r,p);
            
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
    ir = find(reachVec == r,1,'first');
    xr(r) = skm(ir);
    yr(r) = truthAllign.(zField)(ir,3);
end
subplot(2,1,1)
scatter(xr,yr, 100, reaches,'filled')
subplot(2,1,2)
scatter(xr,yr, 100, reaches,'filled')

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position', [0.2, 0.1, 0.6, 0.9])



% %epdf of all node errors
% figure()
% ksdensity(reshape(zErr,[],1))
% hold on
% ksdensity(reshape(z2Err,[],1))
% xlabel('Elevation Error (m)')
% title('EPDF of Node Errors')
% legend('Original Data','Low Rank')

% reach slope errors
simStats.slopeMAE = nanmean(abs(simStats.slopeErr),2);
svdStats.slopeMAE = nanmean(abs(svdStats.slopeErr),2);

simStats.slopeRMSE = sqrt(nanmean(simStats.slopeErr.^2,2));
svdStats.slopeRMSE = sqrt(nanmean(svdStats.slopeErr.^2,2));

simStats.slopeRRMSE = sqrt(nanmean(simStats.relSlopeErr.^2,2)).*100;
svdStats.slopeRRMSE = sqrt(nanmean(svdStats.relSlopeErr.^2,2)).*100;

simStats.slopeRMAE = nanmean(abs(simStats.relSlopeErr),2).*100;
svdStats.slopeRMAE = nanmean(abs(svdStats.relSlopeErr),2).*100;

figure()
subplot(2,2,1)
scatter(simStats.slopeMAE, svdStats.slopeMAE,[],reaches,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated MAE(cm/km)')
ylabel('LRA MAE(cm/km)')
title('MAE')

subplot(2,2,2)
scatter(simStats.slopeRMAE, svdStats.slopeRMAE,[],reaches,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated RMAE(%)')
ylabel('LRA RMAE(%)')
title('RMAE')

subplot(2,2,3)
scatter(simStats.slopeRMSE, svdStats.slopeRMSE,[],reaches,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated RMSE(cm/km)')
ylabel('LRA RMSE(cm/km)')
title('RMSE')

subplot(2,2,4)
% plot(simStats.slopeRRMSE, svdStats.slopeRRMSE,'k.','MarkerSize',20)
scatter(simStats.slopeRRMSE, svdStats.slopeRRMSE,[],reaches,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated RRMSE(%)')
ylabel('LRA RRMSE(%)')
title('RRMSE')

% reach elev. errors
% simStats.reachZMAE = nanmean(abs(simStats.reachZErr),1);
% svdStats.reachZMAE = nanmean(abs(svdStats.reachZErr),1);
% 
% simStats.reachZRMSE = sqrt(nanmean(simStats.reachZErr.^2,1));
% svdStats.reachZRMSE = sqrt(nanmean(svdStats.reachZErr.^2,1));
% 
% figure()
% plot(simStats.reachZRMSE', svdStats.reachZRMSE','k.','MarkerSize',20)
% % plot(abs(simStats.reachZErr), abs(svdStats.reachZErr),'ko')
% scatter1to1(gca);
% title('Reach-level elevation errors')
% xlabel('Simulated (m)')
% ylabel('LRA (m)')

%--------------------------------------------------------------------------


