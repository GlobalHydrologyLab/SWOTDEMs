% DEMErrorEval.m
%
% find std dev of DEM elevations relative to drifter profiles. Reported
% error characteristics are not appropriate for this as they are not
% specific to water surfaces.

clear
close all

load('Sacramento/transformedSacDataV2.mat')

[interpS, interpPts] = unique(drifter(:,4));
interpZ = drifter(interpPts,3);


% Aster
ASTER = nanRows(ASTER);
[~, aster.unqID] = unique(ASTER(:,4));
ASTER = ASTER(aster.unqID,:);
aster.zInterp = interp1(ASTER(:,4),ASTER(:,3),interpS);
aster.bias = nanmean(aster.zInterp - interpZ);
aster.rmBias = aster.zInterp - aster.bias;

aster.resid = aster.rmBias - interpZ;
aster.std = nanstd(aster.resid);


% NED
NED = nanRows(NED);
[~, ned.unqID] = unique(NED(:,4));
NED = NED(ned.unqID,:);
ned.zInterp = interp1(NED(:,4),NED(:,3),interpS);
ned.bias = nanmean(ned.zInterp - interpZ);
ned.rmBias = ned.zInterp - ned.bias;

ned.resid = ned.rmBias - interpZ;
ned.std = nanstd(ned.resid);


% SRTM
SRTM = nanRows(SRTM);
[~, srtm.unqID] = unique(SRTM(:,4));
SRTM = SRTM(srtm.unqID,:);
srtm.zInterp = interp1(SRTM(:,4),SRTM(:,3),interpS);
srtm.bias = nanmean(srtm.zInterp - interpZ);
srtm.rmBias = srtm.zInterp - srtm.bias;

srtm.resid = srtm.rmBias - interpZ;
srtm.std = nanstd(srtm.resid);

%SLM profile
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacSLMProfile.mat')
slm.zInterp = interp1(slm.x,slmProf,interpS);
slm.bias = nanmean(slm.zInterp - interpZ);
slm.rmBias = slm.zInterp - slm.bias;


% c = distinguishable_colors(5);
% c = brewermap(
% set(groot,'defaultAxesColorOrder',c)

skm = interpS./1000;

hold on
plot(skm,aster.rmBias, 'Linewidth',1);
plot(skm,ned.rmBias,'Linewidth',2);
plot(skm,srtm.rmBias,'Linewidth',2);
plot(skm,interpZ,'k','Linewidth',2.5);
plot(skm,slm.rmBias,'Linewidth',2);
legend('ASTER','NED','SRTM','GPS Profile','Simulated SWOT','Location','NorthEast')
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
box on