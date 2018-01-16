% SVD_multiSect_varReach.m
% Extending idea of SVD application to reduce errors in each profile
% individually to work for full river profiles instead of subsections.
% Reaches are defined by equiReach.m, which takes in total number of nodes
% and desired number of nodes per reach.

%   TO DO:
% - look into removing erroneous data from orbit 527 before SVD. The large
%       (and obvious) errors can dominate the decomposition results and
%       make comparisons of errors kind of disingenuous.
% - make sure # of singular values chosen is providing best results.

clear
close all

% % Sac
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
zField = 'geoHeight';
% hard-coded removal of far range data from pass 527
for i = 10:17
    simulated(i).geoHeight(339:487) = NaN;
end

% % Po
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoDataV2.mat')
% zField = 'nHeight';
% % % hard-coded removal of far range data
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
for i = 1:length(simulated)
    reachVec = equiReach(range(nodeRng)+1,nodePerReach);
    simulated(i).reach = reachVec;
    truth(i).reach = reachVec;
end
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
    
    flipFlag=0;
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
    
    %SVD requires nRows>nCol so if this is a problem, flip matrix.
    if size(zresid,2)>size(zresid,1)
        flipFlag = 1;
        zresid = zresid';
    end
    
    [U,S,V] = svd(zresid,0);
    
    %----------------------------------------------------------------------
    % Determine best rank 
    %
    % This section needs to be finished- both break in magnitude and IQR
    % approach perform better than the other in some cases. Maybe some kind
    % of consensus with a third method could work.
    %----------------------------------------------------------------------
    S = diag(S); %extract diagonal values

    % WORKING IDEAS:
    
    %largest decrease in magnitude
%     [~,iSV] = min(diff(S));

    %inter quartile range approach.
    IQR = iqr(S);
    iSV = find(S >= median(S) + 2*IQR,1,'last');
    
    %cumulative percent of sing. val. threshold.
%     p = 0.5;
%     iSV = find(cumsum(S)./sum(S) >= p,1,'first')
    

    S = diag(S); %recreate diagonal matrix
    %----------------------------------------------------------------------


    %now modify S, removing smaller components.
    S2 = S;
    S2(iSV+1:end,:) = 0;
    z2resid = U*S2*V';
    
    %if we had to flip the matrix, flip back.
    if flipFlag
        zresid = zresid';
        z2resid = z2resid';
    end

    z2 = z2resid + polyval(m,s,[],mu);
  
    % join section data for later comparison
    z2All(inSect,~delCol) = z2;
    sAll(inSect,~delCol) = s;
end

skm = nanmean(sAll,2)/1000;
zAll = simAllign.(zField);

%--------------------------------------------------------------------------
% Reach Slopes
%--------------------------------------------------------------------------
reaches = unique(reachVec)';

simSlopeErr = nan(length(reaches),nProf);
SVDSlopeErr = nan(length(reaches),nProf);
reachZErr = nan(length(reaches),nProf);
reachZ2Err = nan(length(reaches),nProf);

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

            simSlopeErr(r,p) = fitSim(1) - fitTruth(1);
            SVDSlopeErr(r,p) = fitSVD(1) - fitTruth(1);
            
            reachAvgZ(r,p) = mean(zAll(inReach,p));
            reachAvgZ2(r,p) = mean(z2All(inReach,p));
            reachAvgTruth(r,p) = mean(truthAllign.(zField)(inReach,p));
            
            reachZErr(r,p) = reachAvgTruth(r,p) - reachAvgZ(r,p);
            reachZ2Err(r,p) = reachAvgTruth(r,p) - reachAvgZ2(r,p);
        end

%         simSlope(r,p) = fitSim(1);
%         SVDSlope(r,p) = fitSVD(1);
%         truthSlope(r,p) = fitTruth(1);
    end

end

%% 

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

%rothko section plot
imagesc(observedBy .* section)
xlabel('Profile')
ylabel('Node')
c = lines;
c = c(1:max(section),:);
colormap([1 1 1; c])

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
% figure()
% ksdensity(reshape(zErr,[],1))
% hold on
% ksdensity(reshape(z2Err,[],1))
% xlabel('Elevation Error (m)')
% title('Empirical PDF of Node Errors')
% legend('Original Data','Low Rank')


%node errors
RMSE = sqrt(nanmean(zErr.^2,1));
RMSESVD = sqrt(nanmean(z2Err.^2,1));
MAE = nanmean(abs(zErr));
MAESVD = nanmean(abs(z2Err));
figure()
bar([RMSE' RMSESVD'],1,'grouped')
% bar([MAE' MAESVD'],1,'grouped')
ylabel('RMSE (m)')
xlabel('Profile Number')
title('Node-level height errors')
legend('Original Data','Low Rank','Location','Northwest')
c = gray;
colormap(c([10,40],:))


% reach slope errors
simSlopeMAE = nanmean(abs(simSlopeErr),1);
SVDSlopeMAE = nanmean(abs(SVDSlopeErr),1);
simSlopeRMSE = sqrt(nanmean(simSlopeErr.^2,1));
SVDSlopeRMSE = sqrt(nanmean(SVDSlopeErr.^2,1));
figure()
% bar([simSlopeMAE' SVDSlopeMAE'] .* 10^5,1,'grouped')
bar([simSlopeRMSE' SVDSlopeRMSE'] .* 10^5,1,'grouped')
ylabel('RMSE (cm/km)')
xlabel('Profile Number')
title('Reach-level slope errors')
legend('Original Data','Low Rank','Location','Northwest')
%     c = lines;
%     colormap(c(1:2,:))
c = gray;
colormap(c([10,40],:))


% reach elev. errors
reachZMAE = nanmean(abs(reachZErr),1);
reachZ2MAE = nanmean(abs(reachZ2Err),1);
reachZRMSE = sqrt(nanmean(reachZErr.^2,1));
reachZ2RMSE = sqrt(nanmean(reachZ2Err.^2,1));
figure()
% bar([simSlopeMAE' SVDSlopeMAE'] .* 10^5,1,'grouped')
bar([reachZRMSE' reachZ2RMSE'] .* 10^5,1,'grouped')
ylabel('RMSE (cm/km)')
xlabel('Profile Number')
title('Reach-level elevation errors')
legend('Original Data','Low Rank','Location','Northwest')
%     c = lines;
%     colormap(c(1:2,:))
c = gray;
colormap(c([10,40],:))

% % R-mode analysis
% R = V*S2;
% figure()
% plot(R(:,1),R(:,2),'r.','MarkerSize',24)

%--------------------------------------------------------------------------


