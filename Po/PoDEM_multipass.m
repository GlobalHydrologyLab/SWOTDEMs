%% PoDEM_multipass.m 
%open SWOT simulator and DEM sampled profiles of the Po River, transform to
%s,n coordinates and compare.
% 
% clear
% 
% %using Legleiter xy2sn.m conversion function
% addpath('cited_functions/riverKrige');
% 
% %load data
% load('/Users/Ted/Documents/MATLAB/Po/PoSimulatorDataV2.mat')
% 
% 
% %% transform
% 
% % transParam = [3 3 31 length(SRTM) 200]'; 
% transParam = [1 3 5 length([truth(1).easting]) 500]'; %for swot sim centerline
% 
% %quick way to change reference cetnerline
% centerline = [truth(1).easting,truth(1).northing];
% 
% 
% %----------------------------------------------------------------------
% %because 0 turns to NaN in the xy2sn function. will not be right if
% %centerline changes.
% % snCoord(1,1)=0;
% %----------------------------------------------------------------------
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
% for i = 1:length(truth)
%     snCoord = xy2sn(centerline, ... 
%         [truth(i).easting, truth(i).northing],transParam);
%     xy2snAscendCheck(snCoord(:,1)); 
%     
%     truth(i).sCoord = snCoord(:,1);
%     truth(i).nCoord = snCoord(:,2);
% end

% % % % % clearvars -except simulated transParam truth clOut
% % % % % save('Po/transformedPoDataV2.mat')

%% Load transformed data from last section
%the coordinate transformation takes far longer than an other process, so
%only run that code when necessary.
% 
clear
% close all

load('Po/transformedPoDataV2.mat')
zField = 'nHeight';


% range=1:3;
% simulated = simulated(range);
% truth = truth(range);

%% trim near and far range
nodeRange = [100;656];
truth = trimFields(truth,nodeRange);
simulated = trimFields(simulated,nodeRange);

%% bias correction
% simAllign = nodeAllign(simulated);
% truthAllign = nodeAllign(truth);
% profMask = ~isnan(simAllign.(zField));
% 
% for i = 1:length(simulated)
%     bias(i) = nanmean(simAllign.(zField)(:,i) - truthAllign.(zField)(:,i));
%     simulated(i).(zField) = simulated(i).(zField) - bias(i);
% end

%% try some clever averaging

[simAvg, profMask] = nodeAvg3_1(simulated, zField);
truthAvg = nodeAvg3_1(truth, zField);%, profMask);


simSmooth = smooth(simAvg.sCoord, simAvg.(zField), 5, 'moving');
truthSmooth = smooth(truthAvg.sCoord, truthAvg.(zField), 5, 'moving');

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


% try slm fit with all data pts and 1/nHeightStd weights. 

%% RMSEs
RMSEslm = sqrt(mean((slmProf - truthAvg.(zField)).^2));
MAEslm = mean(abs(slmProf - truthAvg.(zField)));

RMSEsimAvg = sqrt(mean((simAvg.(zField) - truthAvg.(zField)).^2));
RMSEsmooth = sqrt(nanmean((simSmooth - truthAvg.nHeight).^2));

% for i=1:length(simulated)
%     RMSEsim(i,1) = sqrt(nanmean((simulated(i).(zField) - truth(i).(zField)).^2));
% end


%% plot things
% close all

figure()
plot(truthAvg.sCoord/1000,truthAvg.(zField),'k-','Linewidth',2)
hold on
plot(simAvg.sCoord/1000,simAvg.(zField),'b-','Linewidth',1)
plot(slm.x/1000,slmProf,'r-','LineWidth',2) %slm profile
% plot(slm.knots/1000,knotZ,'b*','MarkerFaceColor','b')
xlabel('Flow distance (km)')
ylabel('Elevation (m)')
hold off

% set(gcf,'Units','normalized','Position', [0.3, 0.2, 0.5, 0.6])

% figure()
% plot(simAvg.sCoord/1000,simAvg.(zField),'r-','Linewidth',2)
% hold on
% plot(truthAvg.sCoord/1000,truthAvg.(zField),'k','Linewidth',2)

