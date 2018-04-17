%% plotProfiles.m
%interpolate DEM derived profiles, compare, plot.

clear
close all

load('smoothingTest/k13_sectMin3.mat')
maxDiff = 0;

%% Po
load('Po/DEMProfiles/dems.mat')
riv = 'Po';

avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
Po.meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(Po.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 750]';
   
[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
xy2snAscendCheck(snCoords(:,1));
Po.skm = snCoords(:,1)./1000;

Po.SRTM = profResample(srtm,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.MERIT = profResample(merit,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.TINITALY = profResample(tinitaly,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.ASTER = profResample(ASTER,clOut,Po.skm,transParam,meanZ,maxDiff);

Po.profNames = fields(Po);
Po.profNames = Po.profNames(3:end);
Po.avgTruth = nanmean(SIMStats.Po.z - SIMStats.Po.zErr,2);
Po.errTable = errorTable(Po,Po.profNames,Po.avgTruth);

%% Sacramento
load('Sacramento/DEMProfiles/dems.mat')
riv = 'Sacramento';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
Sac.meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(Sac.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 200]';

[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
xy2snAscendCheck(snCoords(:,1));
Sac.skm = snCoords(:,1)./1000;

Sac.ASTER = profResample(ASTER,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.NED = profResample(NED,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.SRTM = profResample(SRTM,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.MERIT = profResample(MERIT,clOut,Sac.skm,transParam,meanZ,maxDiff); 
Sac.lidar = profResample(lidar,clOut,Sac.skm,transParam,meanZ,maxDiff); 

%% Tanana
load('Tanana/DEMProfiles/demProfiles.mat')
riv = 'Tanana';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
Tanana.meanProf = SVDStats.(riv).meanProf;
meanZ = nanmean(Tanana.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 200]';

[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
xy2snAscendCheck(snCoords(:,1));
Tanana.skm = snCoords(:,1)./1000;

Tanana.MERIT = profResample(MERIT,clOut,Tanana.skm,transParam,meanZ,maxDiff);
Tanana.ArcticDEM = profResample(ArcticDEM,clOut,Tanana.skm,transParam,meanZ,maxDiff);
Tanana.TanDEMX = profResample(TanDEMX,clOut,Tanana.skm,transParam,meanZ,maxDiff);
Tanana.ASTER = profResample(ASTER,clOut,Tanana.skm,transParam,meanZ,maxDiff);

%% save data

save('Figures/DEMCompare/interpDEMs.mat', 'Po','Sac','Tanana')

%% Plots
clear
close all
load('Figures/DEMCompare/interpDEMs.mat')
LW = 1;
xBuff = 0.025;
yBuff = 0.1;

colors = brewermap(8,'dark2');

figure(1)
set(gcf,'Units','normalized','Position',[0.32461 0.11319 0.31445 0.70556])

subplot(3,1,1)
hold on
plot(Po.skm,Po.SRTM,'Color',colors(4,:),'Linewidth',LW)
plot(Po.skm,Po.MERIT,'Color',colors(6,:),'Linewidth',LW)
plot(Po.skm,Po.ASTER,'Color',colors(3,:),'Linewidth',LW)
plot(Po.skm,Po.TINITALY,'Color',colors(5,:),'Linewidth',LW)
plot(Po.skm,Po.meanProf,'k','Linewidth',2)
centerData(Po.skm,Po.meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('SRTM','MERIT','ASTER','TINITALY','SWOT')
title('Po')
box on


subplot(3,1,2)
hold on
plot(Sac.skm,Sac.SRTM,'Color',colors(4,:),'Linewidth',LW)
plot(Sac.skm,Sac.MERIT,'Color',colors(6,:),'Linewidth',LW)
plot(Sac.skm,Sac.ASTER,'Color',colors(3,:),'Linewidth',LW)
plot(Sac.skm,Sac.NED,'Color',colors(1,:),'Linewidth',LW)
plot(Sac.skm,Sac.lidar,'Color',colors(8,:),'Linewidth',LW)
plot(Sac.skm,Sac.meanProf,'k','Linewidth',2)
centerData(Sac.skm,Sac.meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('SRTM','MERIT','ASTER','NED','LiDAR','SWOT')
title('Sacramento')
box on


subplot(3,1,3)
hold on
plot(Tanana.skm,Tanana.MERIT,'Color',colors(6,:),'Linewidth',LW)
plot(Tanana.skm,Tanana.ASTER,'Color',colors(3,:),'Linewidth',LW)
plot(Tanana.skm,Tanana.ArcticDEM,'Color',colors(2,:),'Linewidth',LW)
plot(Tanana.skm,Tanana.TanDEMX,'Color',colors(7,:),'Linewidth',LW)
plot(Tanana.skm,Tanana.meanProf,'k','Linewidth',2)
centerData(Tanana.skm,Tanana.meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('MERIT','ASTER','ArcticDEM','TanDEM-X','SWOT')
title('Tanana')
box on

%%
function [prof_rs] = profResample(prof,centerline,skm,transParam,meanZ,maxDiff)
%project profile onto SWOT sim centerline, interpolate, remove bias.
prof(:,4:5) = xy2sn(centerline,prof(:,1:2),transParam);  
prof = nanRows(prof(:,1:4)); 
[~,uniqueS] = unique(prof(:,4));
prof_rs = interp1(prof(uniqueS,4)./1000,prof(uniqueS,3),skm);
% prof_rs = slopeConstrain(prof_rs,maxDiff);
prof_rs = prof_rs - (nanmean(prof_rs) - meanZ); 

end

function [t] = errorTable(profiles,profNames,truth)
%structure containing all DEM-derived profiles. compare all with truth and
%compile table.
for i = 1:numel(profNames)
    profErr = profiles.(profNames{i}) - truth;
    MAE(i,1) = nanmean(abs(profErr));
    RMSE(i,1) = sqrt(nanmean(profErr.^2));  
end
t = table(MAE,RMSE,'RowNames',profNames);
end