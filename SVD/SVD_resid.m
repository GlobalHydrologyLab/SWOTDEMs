% SVD_resid.m
% attempting the application of SVD to river data again. This time taking
% after the ideas from the data analysis term paper. try to reduce errors
% in each profile by finding low-rank approximation of all data.

%   TO DO:
% - extend to loop through sections as defined by subsectByObs.m
% -- join those sections
% -- might need to group shorter sections and remove extra observations
% - figure out best number of singular values to use.
% - apply smoothing to each profile

clear
close all

 % Sac
% load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Sacramento/transformedSacDataV2.mat')
% zField = 'geoHeight';
% rng = [2,330]; %upstream
% rng = [475,705]; %downstream

% Po
load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Po/transformedPoData.mat')
zField = 'nHeight';
rng = [100,340];
% 
% Tanana
% load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Tanana/transformedTananaData.mat')
% zField = 'nHeight';
% rng = [503, 624];

simulated = trimFields(simulated,rng);
truth = trimFields(truth,rng);
clearvars -except simulated truth zField

% % hard-coded removal of high mean z days for sacramento
% simulated(16) = [];
% truth(16)=[];
% simulated(6) = [];
% truth(6)=[];

% simAllign = nodeAllign(simulated);
% [sections, observedBy] = subsectByObs(simAllign.(zField));
[simAvg, profMask] = nodeAvg3_1(simulated,zField);
truthAvg = nodeAvg3_1(truth, zField, profMask);

nProf = length(simulated);
polyOrder = 1;

% fit polynomial
s = [simulated.sCoord];
z = [simulated.(zField)];
sv = reshape(s,[],1);
zv = reshape(z,[],1);
[m,~,mu] = polyfit(sv,zv,polyOrder);
zhat = polyval(m,sv,[],mu);
zresid = reshape(zv-zhat, [], nProf);

truthZ = [truth.(zField)];


[U,S,V] = svd(zresid,0);


%--------------------------------------------------------------------------
% This needs evaluated. 
% How to determine best rank- some choices actually increase errors.

S = diag(S); %extract diagonal values

% WORKING IDEAS:

%largest decrease in magnitude
% [~,iSV] = min(diff(S));

%inter quartile range approach.
IQR = iqr(S);
iSV = find(S >= median(S) + 1.5*IQR,1,'last');

S = diag(S); %recreate diagonal matrix
%--------------------------------------------------------------------------


%now modify S, removing smaller components.
S2 = S;
S2(iSV+1:end,:) = 0;
z2 = U*S2*V';

z2 = z2 + polyval(m,s,[],mu);

skm = truthAvg.sCoord/1000;

% smooth?
% for i = 1:nProf
%     z2(:,i) = smooth(z2(:,i),11,'loess');
% end


%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

%original and approx profiles
handle = figure();
subplot(2,1,1);
plot(skm,z)
hold on
% plot(skm,truthAvg.(zField),'k','Linewidth',2)
plot(skm,truth(1).(zField),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Original Simulated Profiles')

subplot(2,1,2);
plot(skm,z2)
hold on
% plot(skm,truthAvg.(zField),'k','Linewidth',2)
plot(skm,truth(1).(zField),'k','Linewidth',2)
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
figure()
ksdensity(reshape(zErr,[],1))
hold on
ksdensity(reshape(z2Err,[],1))
xlabel('Elevation Error (m)')
title('Empirical PDF of Node Errors')
legend('Original Data','Low Rank')

%errors
RMSE = sqrt(mean(zErr.^2,1));
RMSESVD = sqrt(mean(z2Err.^2,1));
MAE = mean(abs(zErr),1);
MAESVD = mean(abs(z2Err),1);
figure()
bar([RMSE' RMSESVD'],1,'grouped')
% bar([MAE' MAESVD'],1,'grouped')
ylabel('RMSE (m)')
xlabel('Profile Number')
title('Comparison of RMSEs')
legend('Original Data','Low Rank','Location','Northwest')
c = lines;
colormap(c(1:2,:))


% % minimum reach length from LeFavour & Alsdorf, 2005.
% Smin = min(abs(diff(truthAvg.(zField))./diff(truthAvg.sCoord)));
% zRL = 2.*std(zErr)./Smin/1000;
% z2RL = 2.*std(z2Err)./Smin/1000;
% if Smin == 0
%     warning('Minimum slope = 0. Can''t calc minimum reach length.')
% else
%     figure()
%     bar([zRL' z2RL'],1,'grouped')
%     ylabel('Minimum Reach Length (km)')
%     xlabel('Profile Number')
%     title('Comparison of Minimum Reach Length')
%     legend('Original Data','Low Rank','Location','Northwest')
%     c = lines;
%     colormap(c(1:2,:))
% end

%R-mode analysis
R = V*S2;
figure()
plot(R(:,1),R(:,2),'r.','MarkerSize',24)




