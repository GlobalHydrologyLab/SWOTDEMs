%% TananaDEM.m 
%open SWOT simulator and DEM sampled profiles of the Tanana, transform to
%s,n coordinates and compare.
% 
% clear
% 
% %using Legleiter xy2sn.m conversion function
% addpath('cited_functions/riverKrige');
% 
% %load data
% load('Tanana/simulated_2cycle.mat')
% load('Tanana/TanDEM-X/TDXProfile.mat')
% load('Tanana/truth.mat')
% 
% 
% %% transform
% 
% transParam = [3 3 11 length(simulated(1).node)*10 500]';
% 
% [snCoord,~,centerlineOut] = xy2sn([TDX(:,1),TDX(:,2)],[TDX(:,1),TDX(:,2)],transParam);
% TDX(:,4:5) = snCoord;
% 
% for i = 1:length(simulated)
%     [snCoord,~,centerlineOut,xyOut] = xy2sn(centerlineOut(:,1:2), ...
%         [simulated(i).easting, simulated(i).northing],transParam);
%     
%     simulated(i).sCoord = snCoord(:,1);
%     simulated(i).ncoord = snCoord(:,2);
%     
%     xy2snAscendCheck(snCoord(:,1));
%     
%     [snCoord,~,centerlineOut,xyOut] = xy2sn(centerlineOut(:,1:2), ...
%         [truth(i).easting, truth(i).northing],transParam);
%     
%     truth(i).sCoord = snCoord(:,1);
%     truth(i).ncoord = snCoord(:,2);
%     
%     xy2snAscendCheck(snCoord(:,1));
% 
% end
% 
% clearvars -except transParam simulated truth TDX
% % % % % save('Tanana/transformedTananaData.mat')

%% Load transformed data from last section
%the coordinate transformation takes far longer than an other process, so
%it makes sense to do it once.

clear
close all

load('Tanana/transformedTananaData.mat')


%% try some clever averaging


simAvg = nodeAvg2(simulated,1);
truthAvg = nodeAvg2(truth,1);
% [truthAvg, wTruth] = nodeAvg2(truth,0.5);
wSim(isnan(wSim)) = 1; %weights = 1/var, so with low num of profiles, var=0 and weights are NaN


simSmooth = smooth(simAvg.sCoord, simAvg.geoHeight, 5, 'moving');
% truthSmooth = smooth(truthAvg.sCoord, truthAvg.geoHeight, 5, 'moving');

%% SLM toolbox 
%trying to automate number of knots and their placement instead of
%informing the number of knots on matching the truth data. Idea is that
%low/high slope areas also have different widths, so placing knots at the
%peaks in widths will approximate the boundaries between low/high slope.

prescription = slmset('Decreasing','on', 'Knots',3, 'Weights',simAvg.nWidth, ...
    'Verbosity',1);
notNaN = ~isnan(simAvg.sCoord); %slm removes these


[~,iKnots] = slidePeaks(simAvg.nWidth(notNaN),101,0.10,0);

locKnots = simAvg.sCoord(iKnots); %s-coordinate location for knots
locKnots(isnan(locKnots)) = []; %delete NaNs
locKnots = [min(simAvg.sCoord); locKnots; max(simAvg.sCoord)];

%high-knot spline
prescription.Knots = locKnots;
slm = slmengine(simAvg.sCoord, simAvg.geoHeight, prescription);
slmProf = slmeval(slm.x,slm);  %evaluate'


% slmTruth = slmengine(truthAvg.sCoord, truthAvg.geoHeight, prescription, 'Weights', wTruth);
% slmProfTruth = slmeval(slmTruth.x,slmTruth);
% 
% % RMSEslm = sqrt(mean((slmProf - truthAvg.geoHeight(notNaN)).^2));
% RMSEslm = sqrt(mean((slmProf - slmProfTruth).^2));

% %% RMSEs
% 
% RMSEsimAvg = sqrt(mean((simAvg.geoHeight - truthAvg.geoHeight).^2));
% RMSEsmooth = sqrt(nanmean((simSmooth - truthSmooth).^2));
% 
% for i=1:length(simulated)
%     RMSE(i,1) = sqrt(nanmean((simulated(i).geoHeight - truth(i).geoHeight).^2));
% end
% 
%  
% %% Slopes
% 
% slopeSLMTruth = diff(slmProfTruth)./diff(slmTruth.x)*-1;
% slopeSLM = diff(slmProf)./diff(slm.x)*-1;
% slopeDiff = slopeSLM-slopeSLMTruth;
% 
%% plot things
close all

% figure(1)
% hold on
% title('pass 249')
% xlabel('flow distance (km)')
% ylabel('geoHeight (m)')
% 
% figure(2)
% hold on
% title('pass 527')
% xlabel('flow distance (km)')
% ylabel('geoHeight (m)')
% 
% for i = 1:length(simulated) 
%     
%     if contains(simulated(i).name, '249')
%         figure(1)
%         plot(simulated(i).sCoord/1000,simulated(i).geoHeight,'r--')
%         plot(truth(i).sCoord/1000,truth(i).geoHeight,'k-')
%        
%     else
%         figure(2)
%         plot(simulated(i).sCoord/1000,simulated(i).geoHeight,'r--')
%         plot(truth(i).sCoord/1000,truth(i).geoHeight,'k-')
% 
%     end   
% end


%profile view
figure(3)
hold on

% plot(ASTER(:,4)/1000,ASTER(:,3),'Linewidth',2) %ASTER
% plot(SRTM(:,4)/1000,SRTM(:,3),'Linewidth',2) %SRTM
% plot(NED(:,4)/1000,NED(:,3),'Linewidth',2) %NED
% plot(drifter(:,4)/1000,drifter(:,3),'.') %drifter

% plot(simAvg.sCoord/1000,simAvg.geoHeight,'LineWidth',2) %averaged output profile
% plot(truthAvg.sCoord/1000,truthAvg.geoHeight,'LineWidth',2) %averaged input profile
plot(TDX(:,4)./1000,TDX(:,3),'LineWidth',2)
plot(slm.x/1000,slmProf,'LineWidth',2) %slm profile



legend('Constrained weighted spline')
title('Median SWOT simulator profiles')
xlabel('Flow distance (km)')
ylabel('Elevation (m)')

% for i = 1:length(boatProfile)
%     plot(boatProfile(i).sCoord/1000,boatProfile(i).height+29)
% end

set(gcf,'Units','normalized','Position', [0.3, 0.2, 0.5, 0.6])

hold off


