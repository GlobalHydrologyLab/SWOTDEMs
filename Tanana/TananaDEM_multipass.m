%% TananaDEM_multipass.m 
%open SWOT simulator and DEM sampled profiles of the Tanana River, transform to
%s,n coordinates and compare.

%% load, transform, and save data
% clear
% 
% %using Legleiter xy2sn.m conversion function
% addpath('cited_functions/riverKrige');
% 
% %load data
% load('Tanana/TananaSimulatorData.mat')
% load('Tanana/InterpTruth.mat')
% %centerline
% load('Tanana/avgTananaCenterline.mat')
% centerline = flip(centerline);
% 
% % transParam = [3 3 31 length(SRTM) 200]'; 
% transParam = [1 3 5 length(centerline(:,1)) 100]'; %for swot sim centerline
% 
% for i = 1:length(simulated)
%     [snCoord,~,clOut] = xy2sn(centerline, ... 
%         [simulated(i).easting, simulated(i).northing],transParam);
%     xy2snAscendCheck(snCoord(:,1)); 
% 
%     simulated(i).sCoord = snCoord(:,1);
%     simulated(i).nCoord = snCoord(:,2);
% end
% 
% %because all truth data are sampled from the same points, we can shortcut
% %and calculate s,n once then assign to all.
% [snCoord,~,clOut] = xy2sn(centerline, ... 
%     [truth(1).easting, truth(1).northing],transParam);
% xy2snAscendCheck(snCoord(:,1));
% 
% for i = 1:length(simulated) 
%     truth(i).sCoord = snCoord(:,1);
%     truth(i).nCoord = snCoord(:,2);
% end
% 
% 
% clearvars -except simulated transParam truth clOut
% save('Tanana/transformedTananaData.mat')

%% Load transformed data from last section
%the coordinate transformation takes far longer than an other process, so
%only run that code when necessary.
% 
clear
close all

load('Tanana/transformedTananaData.mat')
zField = 'nHeight';

%% bias correction

simAllign = nodeAllign(simulated);
truthAllign = nodeAllign(truth);
profMask = ~isnan(simAllign.(zField));

for i = 1:length(simulated)
    bias(i) = nanmean(simAllign.(zField)(:,i) - truthAllign.(zField)(:,i));
    simulated(i).(zField) = simulated(i).(zField) - bias(i);
end
%% try some clever averaging
% nodeRange = [100;407];
% truth = trimFields(truth,nodeRange);
% simulated = trimFields(simulated,nodeRange);

simAvg = nodeAvg3_1(simulated, zField);
truthAvg = nodeAvg3_1(truth, zField, profMask);

simSmooth = smooth(simAvg.sCoord, simAvg.(zField), 10, 'moving');

%% SLM toolbox 
%using the SLM toolbox developed by John D'Errico, posted on MATLAB FEX.
%trying to automate number of knots and their placement instead of
%informing the number of knots on matching the truth data. Idea is that
%low/high slope areas also have different widths, so placing knots at the
%peaks in widths will approximate the low/high slope sequence interval.

prescription = slmset('Decreasing','on','Verbosity',1,'Weights',simAvg.normWeight);
% notNaN = ~isnan(simAvg.sCoord); %slm removes these

% simAvg = trimFields(simAvg, notNaN);
[~,iKnots] = slidePeaks(simAvg.nWidth,0.10,0);

locKnots = simAvg.sCoord(iKnots); %s-coordinate location for knots
locKnots(isnan(locKnots)) = []; %delete NaNs
locKnots = [min(simAvg.sCoord); locKnots; max(simAvg.sCoord)];

prescription.Knots = locKnots;
slm = slmengine(simAvg.sCoord, simAvg.(zField), prescription);
slmProf = slmeval(slm.x,slm);  %evaluate
knotZ = slmeval(slm.knots,slm);

%% RMSEs
RMSEslm = sqrt(nanmean((slmProf - truthAvg.(zField)).^2));
MAEslm = nanmean(abs(slmProf - truthAvg.(zField)));
MAEsim = mean(abs(simAvg.(zField) - truthAvg.(zField)));

RMSEsimAvg = sqrt(nanmean((simAvg.(zField) - truthAvg.(zField)).^2));
RMSEsmooth = sqrt(nanmean((simSmooth - truthAvg.nHeight).^2));

%% Define reaches

nodeRng = [min(simAllign.node(1,:)), max(simAllign.node(end,:))];
nodePerReach = 25;
reachVec = equiReach(range(nodeRng)+1,nodePerReach);

[simulated.reach] = deal(reachVec);
[truth.reach] = deal(reachVec);

%% Reach stats
reaches = unique(reachVec)';

for r = reaches
    inReach = reachVec == r;
    RL(r) = range(simAvg.sCoord(inReach));
    if sum(inReach)>=2
        fitSim = polyfit(simAvg.sCoord(inReach), simAvg.(zField)(inReach),1);
        fitSLM = polyfit(slm.x(inReach),slmProf(inReach),1);
        
        fitTruth = polyfit(truthAvg.sCoord(inReach), truthAvg.(zField)(inReach),1);

        SimSlopeErr(r) = fitTruth(1) - fitSim(1);
        SLMSlopeErr(r) = fitTruth(1) - fitSLM(1);
        SimRelSlopeErr(r) = SimSlopeErr(r) / fitTruth(1) .*100;
        SLMRelSlopeErr(r) = SLMSlopeErr(r) / fitTruth(1) .*100;
    end
end

simRRMSE = sqrt(nanmean(SimRelSlopeErr.^2));
SLMRRMSE = sqrt(nanmean(SLMRelSlopeErr.^2));

%% plot things
close all

figure()
plot(truthAvg.sCoord/1000,truthAvg.(zField),'k-','Linewidth',2)
hold on
plot(slm.x/1000,slmProf,'r-','LineWidth',2) %slm profile
% plot(slm.knots/1000,knotZ,'b*','MarkerFaceColor','b')
xlabel('Flow distance (km)')
ylabel('Elevation (m)')
hold off

%truth inputs, only over extent observed by simulation
figure()
hold on
for i = 1:length(truth)
plot(truth(i).sCoord(profMask(:,i)),truth(i).nHeight(profMask(:,i)),'k-','Linewidth',2)
end

% set(gcf,'Units','normalized','Position', [0.3, 0.2, 0.5, 0.6])

% figure()
% plot(simAvg.sCoord/1000,simAvg.(zField),'r-','Linewidth',2)
% hold on
% plot(truthAvg.sCoord/1000,truthAvg.(zField),'k','Linewidth',2)

