clear
close all

% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/PA.mat')
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/PA_mD0.mat')

nSlopes = 30;
sigmaMultip = 1;
staticCutoff  = 0.3;

trueAvg = slopeConstrain(SIMStats.Sacramento.z - SIMStats.Sacramento.zErr, ... 
    rOpts.Sacramento.maxDiff);
simAvg = SVDStats.Sacramento.meanProf;
s = SIMStats.Sacramento.s;

[Rd,Gd,sId] = relativeSteepness(s.*1000,trueAvg,nSlopes);
trueCutoff = nanmean(Rd) + sigmaMultip*nanstd(Rd);
% trueCutoff = nanmedian(Rd) + sigmaMultip*iqr(Rd);
% trueCutoff = staticCutoff;
steepMaskTrue = Rd>trueCutoff;
steepIdTrue = sId(Rd>trueCutoff);

close all
figure
subplot(2,1,1)
yyaxis left
plot(s,trueAvg,'Linewidth',1)
hold on
plot(s(steepIdTrue),trueAvg(steepIdTrue),'r.','MarkerSize',10)
ylabel('Elevation (m)')

yyaxis right
plot(s(sId),Rd,'Linewidth',1)
hold on
plot(get(gca,'XLim'),[trueCutoff, trueCutoff],'k--')
scaleAxis(gca,1.5)
set(gca,'XTick',0:15:150)
xlabel('Distance (km)')
ylabel('Relative Steepness (m^{-1})')
title('True')

[Rd,Gd,sId] = relativeSteepness(s.*1000,simAvg,nSlopes);
simCutoff = nanmean(Rd) + sigmaMultip*nanstd(Rd);
% simCutoff = nanmedian(Rd) + sigmaMultip*iqr(Rd);
% simCutoff = staticCutoff;
steepMaskSim = Rd>simCutoff;
steepIdSim = sId(steepMaskSim);

subplot(2,1,2)
yyaxis left
plot(s,simAvg,'Linewidth',1)
hold on
plot(s(steepIdSim),simAvg(steepIdSim),'r.','MarkerSize',10)
ylabel('Elevation (m)')

yyaxis right
plot(s(sId),Rd,'Linewidth',1)
hold on
plot(get(gca,'XLim'),[simCutoff, simCutoff],'k--')
scaleAxis(gca,1.5)
set(gca,'XTick',0:15:150)
xlabel('Distance (km)')
ylabel('Relative Steepness (m^{-1})')
title('SWOT sim')
set(gcf,'Position',[555 181 1564 1131])


correct = steepMaskTrue + steepMaskSim == 2;
type1 = steepMaskTrue - steepMaskSim == -1;
type2 = steepMaskTrue - steepMaskSim == 1;


colors = brewermap(3,'*set1') + 0.1; %brighten
colors(2,2) = colors(2,2) - 0.2; %less green

figure
plot(s,trueAvg,'k','Linewidth',2)
hold on
plot(s(sId(correct)),trueAvg(sId(correct)),'.','MarkerSize',15,'Color',colors(1,:))
plot(s(sId(type1)),trueAvg(sId(type1)),'.','MarkerSize',15,'Color',colors(2,:))
plot(s(sId(type2)),trueAvg(sId(type2)),'.','MarkerSize',15,'Color',colors(3,:))

set(gcf,'Position',[359 516 1878 598])
legend('Profile','Correct','Incorrect','Missed')







