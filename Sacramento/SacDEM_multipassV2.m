%% SacDEM.m 
%open SWOT simulator, drifter, and DEM sampled profiles of the Sacramento
%River, transform to s,n coordinates and compare.
% 
% clear
% 
% %using Legleiter xy2sn.m conversion function
% addpath('cited_functions/riverKrige');
% 
% %load data
% load('Sacramento/drawnCenterline.mat') %hand drawn centerline
% load('Sacramento/SWOTSimData/SimTruthV2.mat') %swot sim data
% load('Sacramento/DEMProfiles/SRTM_drawnCenterline.mat') %SRTM
% load('Sacramento/DEMProfiles/ASTER_drawnCenterline.mat') %ASTER
% load('Sacramento/DEMProfiles/NED_drawnCenterline.mat') %NED
% load('Sacramento/SacDrifter.mat') %some drifter data from Tamlin's Frontiers folder
% load('Sacramento/Longitudinalpro/boatProfileStruct.mat') %Longitudinal profiles from Toby
% 
% 
% %% transform
% 
% % transParam = [3 3 31 length(SRTM) 200]'; 
% transParam = [1 3 5 length(SRTM) 200]'; %for swot sim centerline
% 
% %quick way to change reference cetnerline
% drawnCenterline = [truth(1).easting,truth(1).northing];
% 
% %SRTM
% [SRTM(:,4:5),~,centerlineOut] = xy2sn(drawnCenterline(:,1:2),SRTM(:,1:2),transParam);
% xy2snAscendCheck(SRTM(:,4));
% 
% %ASTER
% ASTER(:,4:5) = xy2sn(centerlineOut(:,1:2),ASTER(:,1:2),transParam);
% xy2snAscendCheck(ASTER(:,4));
% 
% %NED
% NED(:,4:5) = xy2sn(centerlineOut(:,1:2),NED(:,1:2),transParam);
% xy2snAscendCheck(NED(:,4));
% 
% %Drifter
% [drifter(:,4:5),~,~,xyOutD,iMissD] = xy2sn(centerlineOut(:,1:2),drifter(:,1:2),transParam);
% xy2snAscendCheck(drifter(:,4));
% 
% 
% %simulator data is processed to the same nodes, so we can do the
% %coordinate conversion once.
% snCoord = xy2sn(centerlineOut(:,1:2), ... 
%         [simulated(1).easting, simulated(1).northing],transParam);
% xy2snAscendCheck(snCoord(:,1)); 

% %----------------------------------------------------------------------
% %because 0 turns to NaN in the xy2sn function. will not be right if
% %centerline changes.
% snCoord(1,1)=0;  
% %----------------------------------------------------------------------
% 
% for i = 1:length(simulated)
%     %assign s,n node coords.
%     simulated(i).sCoord = snCoord(:,1);
%     simulated(i).nCoord = snCoord(:,2);
% 
%     truth(i).sCoord = snCoord(:,1);
%     truth(i).nCoord = snCoord(:,2);
% end
% 
% 
% for i = 1:length(boatProfile)
%     
%     snCoord = xy2sn(centerlineOut(:,1:2), ...
%         [boatProfile(i).easting, boatProfile(i).northing], transParam);
%     
%     boatProfile(i).sCoord = snCoord(:,1);
%     boatProfile(i).nCoord = snCoord(:,2);
%     
% end
% 
% 
% % % % % % % clearvars -except boatProfile drifter simulated SRTM ASTER NED transParam truth
% % % % % % % save('Sacramento/transformedSacDataV2.mat')

%% Load transformed data from last section
%coordinate transformation is very slow 
tic
clear
close all

load('Sacramento/transformedSacDataV2.mat')
zField = 'geoHeight';

%% cleaning up pass 527
% parts of the Sac are outside the swath of pass 526, so there are
% bad/missing values for many fields. Nodes outside the swath can be
% identified by having NaN value for nWidth, so all fields with nWidth =
% NaN are effectively deleted here. Set to NaN for array size continuity.

for i = 1:length(truth)
    %true where widths are missing or height is default -9999 value.
    iBadHeight = isnan(truth(i).nWidth) | truth(i).geoHeight == -9999; 
    
    simulated(i).reach(iBadHeight) = NaN;
    simulated(i).node(iBadHeight) = NaN;
    simulated(i).easting(iBadHeight) = NaN;
    simulated(i).northing(iBadHeight) = NaN;
    simulated(i).nHeight(iBadHeight) = NaN;
    simulated(i).nWidth(iBadHeight) = NaN;
    simulated(i).geoHeight(iBadHeight) = NaN;
    simulated(i).sCoord(iBadHeight) = NaN; 
    simulated(i).nCoord(iBadHeight) = NaN;
    
    
    truth(i).reach(iBadHeight) = NaN;
    truth(i).node(iBadHeight) = NaN;
    truth(i).easting(iBadHeight) = NaN;
    truth(i).northing(iBadHeight) = NaN;
    truth(i).nHeight(iBadHeight) = NaN;
    truth(i).nWidth(iBadHeight) = NaN;
    truth(i).geoHeight(iBadHeight) = NaN;
    truth(i).sCoord(iBadHeight) = NaN; 
    truth(i).nCoord(iBadHeight) = NaN;
end


%% bias correction
simAllign = nodeAllign(simulated);
truthAllign = nodeAllign(truth);
[contSect, profMask] = subsectByObs(simAllign.(zField));

for i = 1:length(simulated)
    bias(i) = nanmean(simAllign.(zField)(:,i) - truthAllign.(zField)(:,i));
    simulated(i).(zField) = simulated(i).(zField) - bias(i);
end


%% ------------------------------------------------------------------------
%% strange problem. 11/3/17
% for some reason, removing the bias results in increased error statistics
% when using nodeAvg2.m but decreased errors with nodeAvg3.m. Need to look
% into the exact differences in these methods, nodeAvg3 was written
% hastily to overcome data organization problems for the Tanana and Po
% river simulations. Should figure out best way to solve 'nodeAvg' and use
% that for all three test cases.
%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------

%% try some clever averaging
% simAvg = nodeAvg2(simulated,0.5);
% truthAvg = nodeAvg2(truth,0.5);
simAvg = nodeAvg3_1(simulated,'geoHeight');
truthAvg = nodeAvg3_1(truth,'geoHeight');

simSmooth = smooth(simAvg.sCoord, simAvg.geoHeight, 5, 'moving');
truthSmooth = smooth(truthAvg.sCoord, truthAvg.geoHeight, 5, 'moving');

%% SLM toolbox 
%using the SLM toolbox developed by John D'Errico, posted on MATLAB FEX.
%trying to automate number of knots and their placement instead of
%informing the number of knots on matching the truth data. Idea is that
%low/high slope areas also have different widths, so placing knots at the
%peaks in widths will approximate the low/high slope sequence interval.

prescription = slmset('Decreasing','on','Verbosity',1,'Weights',simAvg.normWeight);

[~,iKnots] = slidePeaks(simAvg.nWidth,0.10,0);

locKnots = simAvg.sCoord(iKnots); %s-coordinate location for knots
locKnots(isnan(locKnots)) = []; %delete NaNs
locKnots = [min(simAvg.sCoord); locKnots; max(simAvg.sCoord)];

prescription.Knots = locKnots;
slm = slmengine(simAvg.sCoord, simAvg.geoHeight, prescription);
slmProf = slmeval(slm.x,slm);  %evaluate
knotZ = slmeval(slm.knots,slm);

%% RMSEs
RMSEslm = sqrt(mean((slmProf - truthAvg.geoHeight).^2));
MAEslm = mean(abs(slmProf - truthAvg.geoHeight));

RMSEsimAvg = sqrt(mean((simAvg.geoHeight - truthAvg.geoHeight).^2));
RMSEsmooth = sqrt(nanmean((simSmooth - truthSmooth).^2));

for i=1:length(simulated)
    RMSEsim(i,1) = sqrt(nanmean((simulated(i).geoHeight - truth(i).geoHeight).^2));
    RMSEavg(i,1) = sqrt(nanmean((simAvg.geoHeight - truth(i).geoHeight).^2));
end

 
%% Slopes

% slopeSLMTruth = diff(slmProfTruth)./diff(slmTruth.x)*-1;
slopeTruth = diff(truthAvg.sCoord)./diff(truthAvg.geoHeight)*-1;
slopeSLM = diff(slmProf)./diff(slm.x)*-1;
slopeDiff = slopeSLM-slopeTruth;

%% plot things
close all
% 
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
% figure(1)
% h = findobj(gca,'Type','line');
% legend(h(1:2),'model','simulated');
% 
% figure(2)
% h = findobj(gca,'Type','line');
% legend(h(1:2),'model','simulated');



%profile view
figure()
hold on

% plot(ASTER(:,4)/1000,ASTER(:,3),'Linewidth',2)
% plot(SRTM(:,4)/1000,SRTM(:,3),'Linewidth',2)
% plot(NED(:,4)/1000,NED(:,3),'Linewidth',2) 
% plot(drifter(:,4)/1000,drifter(:,3)+29,'Linewidth',2)
% xlabel('flow distance (km)')
% ylabel('elevation (m)')
% box on
% legend('ASTER','SRTM','NED','Drifter')
% pbaspect([4 3 1])


% 
% plot(simAvg.sCoord/1000,simAvg.geoHeight,'LineWidth',2) %averaged output profile
% plot(truthAvg.sCoord/1000,truthAvg.geoHeight,'LineWidth',2) %averaged input profile
% plot(truth(1).sCoord/1000,[truth.geoHeight])
plot(truthAvg.sCoord/1000,truthAvg.geoHeight,'k','LineWidth',2) %averaged input profile
plot(slm.x/1000,slmProf,'r-','LineWidth',2) %slm profile
plot(simAvg.sCoord/1000,simAvg.geoHeight,'-')

legend('Simulator input node median','Constrained weighted spline')
title('Median SWOT simulator profiles')
xlabel('Flow distance (km)')
ylabel('Elevation (m)')

% for i = 1:length(boatProfile)
%     plot(boatProfile(i).sCoord/1000,boatProfile(i).height+29)
% end

% set(gcf,'Units','normalized','Position', [0.3, 0.2, 0.5, 0.6])

hold off

toc
