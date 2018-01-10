% SVD_multiReach.m
% Extending idea of SVD application to reduce errors in each profile
% individually to work for full river profiles instead of subsections. For
% now, using reaches as defined by Renato.

%   TO DO:
% - look into removing erroneous data from orbit 527 before SVD. The large
%       (and obvious) errors can dominate the decomposition results and
%       make comparisons of errors kind of disingenuous.
% - make sure # of singular values chosen is providing best results.
% - compare results to a svds() methods with reaches instead of sections.

clear
close all

%Sac
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
zField = 'geoHeight';

% hard-coded removal of far range data from pass 527
for i = 10:17
    simulated(i).geoHeight(339:487) = NaN;
end
% % hard-coded removal of high mean z days for sacramento
simulated(16) = [];
truth(16)=[];
simulated(6) = [];
truth(6)=[];


clearvars -except simulated truth zField



%Group data into matrices without gaps, intelligently deleting data so that
%all sections have >= sectMin rows
sectMin = 25;
simAllign = nodeAllign(simulated);
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
    
    simReach = trimFields(simulated,inSect);
    truthReach = trimFields(truth,inSect);
    
    %assemble full rectangular matrices of s,z data 
    z = zArray(inSect,:);
    [z, delCol] = nanRows(z,2);
    s = simAllign.sCoord(inSect,~delCol);
    
    
    % fit polynomial to all data in reach
    nProf = size(z,2);
    polyOrder = 1;
    sv = reshape(s,[],1);
    zv = reshape(z,[],1);
    [m,~,mu] = polyfit(sv,zv,polyOrder);
    zhat = polyval(m,s,[],mu);
    zresid = reshape(z-zhat, [], nProf);
    
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
%     [~,iSV] = min(diff(S))

    %inter quartile range approach.
    IQR = iqr(S);
    iSV = find(S >= median(S) + 2*IQR,1,'last');
    
    %cumulative percent of sing. val. threshold.
%     p = 0.4;
%     iSV = find(cumsum(S)./sum(S) >= p,1,'first')
    

    S = diag(S); %recreate diagonal matrix
    %----------------------------------------------------------------------


    %now modify S, removing smaller components.
    S2 = S;
    S2(iSV+1:end,:) = 0;
    z2resid = U*S2*V';

    z2 = z2resid + polyval(m,s,[],mu);
  
    % join section data for later comparison
    z2All(inSect,~delCol) = z2;
    sAll(inSect,~delCol) = s;
end

skm = nanmean(sAll,2)/1000;
zAll = simAllign.(zField);
truthZ = truthAllign.(zField);

%--------------------------------------------------------------------------
% Reach Slopes
%--------------------------------------------------------------------------
hasReaches = isfield(simulated,'reach');

if hasReaches
    reaches = unique(simulated(1).reach)';

    for r = reaches
        for p = 1:nProf
            inReach = simAllign.reach(:,p) == r;
            fitSim = polyfit(sAll(inReach,p),zAll(inReach,p),1);
            fitSVD = polyfit(sAll(inReach,p),z2All(inReach,p),1);

            %calc reach length
            RL(r,1) = range(sAll(inReach,1));
            
            inReach = truth(p).reach == r;
            fitTruth = polyfit(truth(p).sCoord(inReach), ... 
                truth(p).(zField)(inReach),1);

            simSlopeErr(r,p) = fitSim(1) - fitTruth(1);
            SVDSlopeErr(r,p) = fitSVD(1) - fitTruth(1);

    %         simSlope(r,p) = fitSim(1);
    %         SVDSlope(r,p) = fitSVD(1);
    %         truthSlope(r,p) = fitTruth(1);
        end


    end
end
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


%original and approx profiles
handle = figure();
subplot(2,1,1);
plot(skm,zAll)
hold on
plot(skm,truth(3).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Original Simulated Profiles')

for r = reaches
    ir = find(simulated(1).reach == r,1,'last');
    xr = simulated(1).sCoord(ir);
    yr = truth(3).(zField)(ir);
    plot([xr xr]./1000,[yr yr], 'r*','LineWidth',2)
end

subplot(2,1,2);
plot(skm,z2All)
hold on
plot(skm,truth(3).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Low-Rank Approximation Profiles')

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position', [0.2, 0.1, 0.6, 0.9])

for r = reaches
    ir = find(simulated(1).reach == r,1,'last');
    xr = simulated(1).sCoord(ir);
    yr = truth(3).(zField)(ir);
    plot([xr xr]./1000,[yr yr], 'r*','LineWidth',2)
end

%singular values
figure()
bar(diag(S))
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
title('Singular Values of Elevation Residuals')


zErr = zAll - truthZ;
z2Err = z2All - truthZ;

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
MAE = nanmean(abs(zErr),1);
MAESVD = nanmean(abs(z2Err),1);
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
%     bar([simSlopeMAE' SVDSlopeMAE'] .* 10^5,1,'grouped')
bar([simSlopeRMSE' SVDSlopeRMSE'] .* 10^5,1,'grouped')
ylabel('RMSE (cm/km)')
xlabel('Profile Number')
title('Reach-level slope errors')
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


