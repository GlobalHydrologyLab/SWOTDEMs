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
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
zField = 'nWidth';
% rng = [2,330]; %upstream
rng = [475,705]; %downstream

% % Po
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPoData.mat')
% zField = 'nHeight';
% rng = [100,340];
% 
% Tanana
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/transformedTananaData.mat')
% zField = 'nHeight';
% rng = [140,458];

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
% Determine best rank 
%
% This section needs to be finished- magnitude option increases errors in
% some cases. 
%--------------------------------------------------------------------------

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
z2resid = U*S2*V';

z2 = z2resid + polyval(m,s,[],mu);

skm = truthAvg.sCoord/1000;


%--------------------------------------------------------------------------
% Compare reach slopes
%--------------------------------------------------------------------------
hasReaches = isfield(simulated,'reach');

if hasReaches
    reaches = unique(simulated(1).reach)';

    for r = reaches

        for p = 1:nProf
            inReach = simulated(p).reach == r;
            fitSim = polyfit(s(inReach,p),z(inReach,p),1);
            fitSVD = polyfit(s(inReach,p),z2(inReach,p),1);

            inReach = truth(p).reach == r;
            fitTruth = polyfit(truth(p).sCoord(inReach), ... 
                truth(p).(zField)(inReach),1);

            simSlopeErr(r,p) = fitSim(1) - fitTruth(1);
            SVDSlopeErr(r,p) = fitSVD(1) - fitTruth(1);

    %         simSlope(r,p) = fitSim(1);
    %         SVDSlope(r,p) = fitSVD(1);
    %         truthSlope(r,p) = fitTruth(1);
        end

        %calc reach length
        RL(r,1) = range(s(inReach,1));

    end
end
%--------------------------------------------------------------------------


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


%node errors
RMSE = sqrt(mean(zErr.^2,1));
RMSESVD = sqrt(mean(z2Err.^2,1));
figure()
bar([RMSE' RMSESVD'],1,'grouped')
ylabel('RMSE (m)')
xlabel('Profile Number')
title('Node-level height errors')
legend('Original Data','Low Rank','Location','Northwest')
c = lines;
colormap(c(1:2,:))


%reach slope errors
if hasReaches
%     simSlopeMAE = mean(abs(simSlopeErr),1);
%     SVDSlopeMAE = mean(abs(SVDSlopeErr),1);
    simSlopeRMSE = sqrt(mean(simSlopeErr.^2,1));
    SVDSlopeRMSE = sqrt(mean(SVDSlopeErr.^2,1));
    figure()
%     bar([simSlopeMAE' SVDSlopeMAE'] .* 10^5,1,'grouped')
    bar([simSlopeRMSE' SVDSlopeRMSE'] .* 10^5,1,'grouped')
    ylabel('RMSE (cm/km)')
    xlabel('Profile Number')
    title('Reach-level slope errors')
    legend('Original Data','Low Rank','Location','Northwest')
    c = lines;
    colormap(c(1:2,:))
end

% % R-mode analysis

% R = V*S2;
% figure()
% plot(R(:,1),R(:,2),'r.','MarkerSize',24)

%--------------------------------------------------------------------------


