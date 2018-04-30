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
% Po.meanProf = SVDStats.(riv).meanProf;
Po.meanProf = slopeConstrain(SIMStats.(riv).z, 0.01);
Po.avgTruth = nanmean(SIMStats.Po.z - SIMStats.Po.zErr,2);
Po.SWOT = SVDStats.(riv).meanProf;
meanZ = nanmean(Po.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 750]';
   
[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
xy2snAscendCheck(snCoords(:,1));
Po.skm = snCoords(:,1)./1000;
Po = clean_skm(Po);

Po.SRTM = profResample(srtm,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.MERIT = profResample(merit,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.ASTER = profResample(ASTER,clOut,Po.skm,transParam,meanZ,maxDiff);
Po.TINITALY = profResample(tinitaly,clOut,Po.skm,transParam,meanZ,maxDiff);

Po.DEMNames = fields(Po);
Po.DEMNames = Po.DEMNames([3,5:end]);
Po.errTable = errorTable(Po,Po.DEMNames,Po.avgTruth,Po.skm);

%% Sacramento
load('Sacramento/DEMProfiles/dems.mat')
riv = 'Sacramento';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
% Sac.meanProf = SVDStats.(riv).meanProf;
Sac.meanProf = slopeConstrain(SIMStats.(riv).z, 0.01);
Sac.avgTruth = nanmean(SIMStats.(riv).z - SIMStats.(riv).zErr,2);
Sac.SWOT = SVDStats.(riv).meanProf;
meanZ = nanmean(Sac.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 200]';

[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
xy2snAscendCheck(snCoords(:,1));
Sac.skm = snCoords(:,1)./1000;
Sac = clean_skm(Sac);

Sac.SRTM = profResample(SRTM,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.MERIT = profResample(MERIT,clOut,Sac.skm,transParam,meanZ,maxDiff); 
Sac.ASTER = profResample(ASTER,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.NED = profResample(NED,clOut,Sac.skm,transParam,meanZ,maxDiff);
Sac.lidar = profResample(lidar,clOut,Sac.skm,transParam,meanZ,maxDiff);

Sac.DEMNames = fields(Sac);
Sac.DEMNames = Sac.DEMNames([3,5:end]);
Sac.errTable = errorTable(Sac,Sac.DEMNames,Sac.avgTruth,Sac.skm);

%% Tanana
load('Tanana/DEMProfiles/demProfiles.mat')
riv = 'Tanana';
avgCenter = [nanmean(SVDStats.(riv).x,2),nanmean(SVDStats.(riv).y,2)];
% Tan.meanProf = SVDStats.(riv).meanProf;
Tan.meanProf = slopeConstrain(SIMStats.(riv).z, 0.005=);
Tan.avgTruth = nanmean(SIMStats.(riv).z - SIMStats.(riv).zErr,2);
Tan.SWOT = SVDStats.(riv).meanProf;
meanZ = nanmean(Tan.meanProf); %for bias removal
transParam = [1 3 5 length(avgCenter)*10 200]';

[snCoords,~,clOut] = xy2sn(avgCenter,avgCenter,transParam);
Tan.skm = snCoords(:,1)./1000;
Tan = clean_skm(Tan);

Tan.MERIT = profResample(MERIT,clOut,Tan.skm,transParam,meanZ,maxDiff);
Tan.ASTER = profResample(ASTER,clOut,Tan.skm,transParam,meanZ,maxDiff);
Tan.ArcticDEM = profResample(ArcticDEM,clOut,Tan.skm,transParam,meanZ,maxDiff);
Tan.TanDEMX = profResample(TanDEMX,clOut,Tan.skm,transParam,meanZ,maxDiff);

Tan.DEMNames = fields(Tan);
Tan.DEMNames = Tan.DEMNames([3,5:end]);

%prepare truth data
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/gpsProfiles/gpsProf.mat')
gps(:,4:5) = xy2sn(clOut,[gps(:,1),gps(:,2)],transParam);
gps = nanRows(gps);
[~,gpsUnique] = unique(gps(:,4));
gps = gps(gpsUnique,:);
gps(:,4) = gps(:,4)./1000; %km 
Tan.errTable = errorTable(Tan,Tan.DEMNames,gps(:,3),gps(:,4));

%% save data

% save('Figures/DEMCompare/interpConstrainDEMs.mat', 'Po','Sac','Tan')

%% Plots
% clear
close all
% load('Figures/DEMCompare/interpConstrainDEMs.mat')
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
plot(Tan.skm,Tan.MERIT,'Color',colors(6,:),'Linewidth',LW)
plot(Tan.skm,Tan.ASTER,'Color',colors(3,:),'Linewidth',LW)
plot(Tan.skm,Tan.ArcticDEM,'Color',colors(2,:),'Linewidth',LW)
plot(Tan.skm,Tan.TanDEMX,'Color',colors(7,:),'Linewidth',LW)
plot(Tan.skm,Tan.meanProf,'k','Linewidth',2)
centerData(Tan.skm,Tan.meanProf,xBuff,yBuff);
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
legend('MERIT','ASTER','ArcticDEM','TanDEM-X','SWOT')
title('Tanana')
box on

%%
function[prof] = clean_skm(prof)
iNaN = find(isnan(prof.skm));
prof.meanProf(iNaN) = [];
prof.skm(iNaN) = [];
prof.avgTruth(iNaN) = [];
prof.SWOT(iNaN) = [];
end


function [prof_rs] = profResample(prof,centerline,skm,transParam,meanZ,maxDiff)
%project profile onto SWOT sim centerline, interpolate, remove bias.
prof(:,4:5) = xy2sn(centerline,prof(:,1:2),transParam); 
prof = nanRows(prof(:,1:4)); 

ds = abs(nanmean(diff(prof(:,4))));
window = 200/ds;
window = ceil(window); %ensure >1
if ~mod(window,2) %ensure odd.
    window = window - 1;
end
if window ~= 1
    prof(:,3) = smooth(prof(:,3),window);
end

[~,uniqueS] = unique(prof(:,4));
prof_rs = interp1(prof(uniqueS,4)./1000,prof(uniqueS,3),skm);
prof_rs = slopeConstrain(prof_rs,maxDiff);
prof_rs = prof_rs - (nanmean(prof_rs) - meanZ); 

% prof = prof(uniqueS,:);
% prof(:,3) = slopeConstrain(prof(:,3),maxDiff);
% prof(:,3) = prof(:,3) - (nanmean(prof(:,3)) - meanZ);
% 
% prof_rs = prof;
end

function [t] = errorTable(profiles,DEMNames,truth,truthS)
%structure containing all DEM-derived profiles. compare all with truth and
%compile table.

for i = 1:numel(DEMNames)
    truth_rs = interp1(truthS,truth,profiles.skm);
    profErr = profiles.(DEMNames{i}) - truth_rs;
    profErr = profErr - nanmean(profErr(:));
    MAE(i,1) = nanmean(abs(profErr));
    RMSE(i,1) = sqrt(nanmean(profErr.^2));  
end
t = table(MAE,RMSE,'RowNames',DEMNames);
end