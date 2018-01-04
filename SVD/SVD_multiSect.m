% SVD_multiReach.m
% Extending idea of SVD application to reduce errors in each profile
% individually to work for full river profiles instead of subsections. For
% now, using reaches as defined by Renato.

%   TO DO:
% - extend to loop through reaches
% -- join reaches, not necessary but provides nice visual.
% - make sure # of singular values chosen is providing best results.

clear
close all

 % Sac
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
zField = 'geoHeight';
clearvars -except simulated truth zField

%sectify
simAllign = nodeAllign(simulated);
simAllign.(zField)(simAllign.(zField)==-9999) = NaN;
truthAllign = nodeAllign(truth);

sectMin = 25;
[section, zArray] = subsectByObs(simAllign.(zField),sectMin);
observedBy = ~isnan(zArray);

%rothko plot
imagesc(observedBy .* section)
xlabel('Profile')
ylabel('Node')
c = lines;
c = c(1:max(section),:);
colormap([1 1 1; c])


for r = min(section):2%max(section)
    
    inSect = section == r;
    
    simReach = trimFields(simulated,inSect);
    truthReach = trimFields(truth,inSect);

%     [simAvg, profMask] = nodeAvg3_1(simReach,zField);
%     truthAvg = nodeAvg3_1(truthReach, zField, profMask);
    
    %assemble full rectangular matricies of s,z data 
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
    % This section needs to be finished- magnitude option increases errors
    % in some cases.
    %----------------------------------------------------------------------

    S = diag(S); %extract diagonal values

    % WORKING IDEAS:

    %largest decrease in magnitude
    % [~,iSV] = min(diff(S));

    %inter quartile range approach.
    IQR = iqr(S);
    iSV = find(S >= median(S) + 1.5*IQR,1,'last');
%     iSV=2;

    S = diag(S); %recreate diagonal matrix
    %----------------------------------------------------------------------


    %now modify S, removing smaller components.
    S2 = S;
    S2(iSV+1:end,:) = 0;
    z2resid = U*S2*V';

    z2 = z2resid + polyval(m,s,[],mu);

%     skm = truthAvg.sCoord/1000;
%     truthZ = [truthReach.(zField)];

end

%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

%original and approx profiles
handle = figure();
subplot(2,1,1);
plot(skm,z)
hold on
% plot(skm,truthAvg.(zField),'k','Linewidth',2)
plot(skm,truthReach(1).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Original Simulated Profiles')

subplot(2,1,2);
plot(skm,z2)
hold on
% plot(skm,truthAvg.(zField),'k','Linewidth',2)
plot(skm,truthReach(1).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Low-Rank Approximation Profiles')

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position', [0.2, 0.1, 0.6, 0.9])


%singular values
figure()
bar(diag(S))
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
title('Singular Values of Elevation Residuals')
% set(gca,'YScale','log')

zErr = z - truthZ;
z2Err = z2 -truthZ;

%epdf of all node errors
% figure()
% ksdensity(reshape(zErr,[],1))
% hold on
% ksdensity(reshape(z2Err,[],1))
% xlabel('Elevation Error (m)')
% title('Empirical PDF of Node Errors')
% legend('Original Data','Low Rank')


%node errors
RMSE = sqrt(mean(zErr.^2,1));
RMSESVD = sqrt(mean(z2Err.^2,1));
figure()
bar([RMSE' RMSESVD'],1,'grouped')
ylabel('RMSE (m)')
xlabel('Profile Number')
title('Node-level height errors')
legend('Original Data','Low Rank','Location','Northwest')
c = gray;
colormap(c([10,40],:))


%reach slope errors
%     simSlopeMAE = mean(abs(simSlopeErr),1);
%     SVDSlopeMAE = mean(abs(SVDSlopeErr),1);
%     simSlopeRMSE = sqrt(mean(simSlopeErr.^2,1));
%     SVDSlopeRMSE = sqrt(mean(SVDSlopeErr.^2,1));
%     figure()
% %     bar([simSlopeMAE' SVDSlopeMAE'] .* 10^5,1,'grouped')
%     bar([simSlopeRMSE' SVDSlopeRMSE'] .* 10^5,1,'grouped')
%     ylabel('RMSE (cm/km)')
%     xlabel('Profile Number')
%     title('Reach-level slope errors')
%     legend('Original Data','Low Rank','Location','Northwest')
% %     c = lines;
% %     colormap(c(1:2,:))
%     c = gray;
%     colormap(c([10,40],:))

% % R-mode analysis

% R = V*S2;
% figure()
% plot(R(:,1),R(:,2),'r.','MarkerSize',24)

%--------------------------------------------------------------------------


