%% Load transformed data
%the coordinate transformation takes far longer than an other process, so
%only run that code when necessary.
% 
clear
close all

load('Po/transformedPoData.mat')
zField = 'nHeight';

tic

nTotProf = length(simulated);
RMSEslmSUM = zeros(nTotProf,1);
MAEslmSUM = zeros(nTotProf,1);
RMSEsimAvgSUM = zeros(nTotProf,1);
nCombRec = zeros(nTotProf,1);

allSimulated = simulated;
allTruth = truth;

minProf = 1;
maxProf = 2%nTotProf;

for numProfiles = minProf:maxProf
    disp(numProfiles)
    maxComb = 10;
    
    %generate all possible combinations of the profiles for current
    %numProfiles.
    allComb = combnk((1:nTotProf),numProfiles);
    nComb = size(allComb,1);
    
    %randomly select maxComb combinations from all valid combinations to
    %reduce run time.
    if nComb > maxComb
        rId = randi(nComb,maxComb,1);
        allComb = allComb(rId,:);
    end
    nComb = size(allComb,1);
    nCombRec(numProfiles) = nComb;
        
    for j = 1:nComb
        
        simulated = allSimulated(allComb(j,:));
        truth = allTruth(allComb(j,:));


        %% try some clever averaging
        nodeRange = [100;407];
        truth = trimFields(truth,nodeRange);
        simulated = trimFields(simulated,nodeRange);

        [simAvg, profMask] = nodeAvg3_1(simulated, zField);
        truthAvg = nodeAvg3_1(truth, zField, profMask);

        simSmooth = smooth(simAvg.sCoord, simAvg.(zField), 5, 'moving');
        truthSmooth = smooth(truthAvg.sCoord, truthAvg.(zField), 5, 'moving');

        %% SLM toolbox 
        %using the SLM toolbox developed by John D'Errico, posted on MATLAB FEX.
        %trying to automate number of knots and their placement instead of
        %informing the number of knots on matching the truth data. Idea is that
        %low/high slope areas also have different widths, so placing knots at the
        %peaks in widths will approximate the low/high slope sequence interval.

        prescription = slmset('Decreasing','on','Verbosity',0,'Weights',simAvg.normWeight);
        notNaN = ~isnan(simAvg.sCoord); %slm removes these

        simAvg = trimFields(simAvg, notNaN);
        [~,iKnots] = slidePeaks(simAvg.nWidth,0.10,0);

        locKnots = simAvg.sCoord(iKnots); %s-coordinate location for knots
        locKnots(isnan(locKnots)) = []; %delete NaNs
        locKnots = [min(simAvg.sCoord); locKnots; max(simAvg.sCoord)];

        prescription.Knots = locKnots;
        slm = slmengine(simAvg.sCoord, simAvg.(zField), prescription);
        slmProf = slmeval(slm.x,slm);  %evaluate
        knotZ = slmeval(slm.knots,slm);
    
        %% RMSEs
        
        RMSEslmSUM(numProfiles) = RMSEslmSUM(numProfiles) + sqrt(nanmean((slmProf - truthAvg.(zField)(notNaN)).^2));
        MAEslmSUM(numProfiles,1) = MAEslmSUM(numProfiles,1) + nanmean(abs(slmProf - truthAvg.(zField)(notNaN)));

        RMSEsimAvgSUM(numProfiles) =  RMSEsimAvgSUM(numProfiles) +  ...
            sqrt(nanmean((simAvg.(zField) - truthAvg.(zField)).^2));
        
    end
    RMSEslm(numProfiles) = RMSEslmSUM(numProfiles) / nComb;
    MAEslm(numProfiles,1) = MAEslmSUM(numProfiles,1) /nComb;
    RMSEsimAvg(numProfiles) =  RMSEsimAvgSUM(numProfiles) /nComb;
    
end

RMSEslm(1:minProf-1)=NaN;
MAEslm(1:minProf-1)=NaN;
RMSEsimAvg(1:minProf-1)=NaN;


toc
%% plots 
close all
plot(RMSEslm,'k','Linewidth',2)
hold on
ylabel('RMSE (m)')
xlabel('# of SWOT observations')
title('Average profile errors')
hold off
% 
% 
% figure()
% histogram(simAvg.geoHeight(notNaN) - truthAvg.geoHeight(notNaN))
% hold on
% histogram(slmProf - truthAvg.geoHeight(notNaN))

