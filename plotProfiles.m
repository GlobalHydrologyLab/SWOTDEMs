%% plotProfiles.m
%interpolate DEM derived profiles, compare, plot.

clear
close all


% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/k_13.mat')
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/k13_sectMin3.mat')
maxDiff = 0;
demLW = 1;
xBuff = 0.025;
yBuff = 0.1;


colors = brewermap(8,'dark2');
set(gca,'ColorOrder',colors)
set(groot,'DefaultAxesColorOrder',colors);

figure
set(gcf,'Units','normalized','Position',[0.32461 0.11319 0.31445 0.70556])


%% Po
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/DEMProfiles/dems.mat')
riv = 'Po';

avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
skm = cumDist(avgCenter(:,1),avgCenter(:,2))./1000;
meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter) 500]';
   

srtm_rs = profResample(srtm,avgCenter,skm,transParam,meanZ);
merit_rs = profResample(merit,avgCenter,skm,transParam,meanZ);
tinitaly_rs = profResample(tinitaly,avgCenter,skm,transParam,meanZ);
%ASTER


subplot(3,1,1)
plot(skm,srtm_rs,'Linewidth',demLW)
hold on
plot(skm,merit_rs,'Linewidth',demLW)
plot(skm,tinitaly_rs,'Linewidth',demLW)
plot(skm,meanProf,'k','Linewidth',2)

centerData(skm,meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('SRTM','MERIT','TINITALY','SWOT')
title('Po')

%% Sacramento
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/DEMProfiles/dems.mat')
riv = 'Sacramento';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
skm = cumDist(avgCenter(:,1),avgCenter(:,2))./1000;
meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter) 500]';


ASTER_rs = profResample(ASTER,avgCenter,skm,transParam,meanZ);
NED_rs = profResample(NED,avgCenter,skm,transParam,meanZ);
SRTM_rs = profResample(SRTM,avgCenter,skm,transParam,meanZ); 
lidar_rs = profResample(lidar,avgCenter,skm,transParam,meanZ); 
%MERIT

subplot(3,1,2)
plot(skm,ASTER_rs,'Linewidth',0.5)
hold on
plot(skm,NED_rs,'Linewidth',demLW)
plot(skm,SRTM_rs,'Linewidth',demLW)
plot(skm,lidar_rs,'Linewidth',demLW)
plot(skm,meanProf,'k','Linewidth',2)

centerData(skm,meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('ASTER','NED','SRTM','LiDAR','SWOT')
title('Sacramento')

%% Tanana
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/DEMProfiles/demProfiles.mat')
riv = 'Tanana';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
skm = cumDist(avgCenter(:,1),avgCenter(:,2))./1000;
meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter) 500]';


MERIT_rs = profResample(MERIT,avgCenter,skm,transParam,meanZ);
ArcticDEM_rs = profResample(ArcticDEM,avgCenter,skm,transParam,meanZ);
TanDEMX_rs = profResample(TanDEMX,avgCenter,skm,transParam,meanZ);

subplot(3,1,3)
plot(skm,MERIT_rs,'Linewidth',demLW)
hold on
plot(skm,ArcticDEM_rs,'Linewidth',demLW)
plot(skm,TanDEMX_rs,'Linewidth',demLW)
plot(skm,meanProf,'k','Linewidth',2)

centerData(skm,meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('MERIT','ArcticDEM','TanDEM-X','SWOT')
title('Tanana')


%%

function [prof_rs] = profResample(prof,centerline,skm,transParam,meanZ)
%project profile onto SWOT sim centerline, interpolate, remove bias.

prof(:,4:5) = xy2sn(centerline,prof(:,1:2),transParam);  
if isnan(prof(1,4)); prof(1,4) = 0; end
prof = nanRows(prof);
[~,uniqueS] = unique(prof(:,4));
prof_rs = interp1(prof(uniqueS,4)./1000,prof(uniqueS,3),skm);
% prof_rs = slopeConstrain(prof_rs,maxDiff);
prof_rs = prof_rs - (nanmean(prof_rs) - meanZ); 

end
