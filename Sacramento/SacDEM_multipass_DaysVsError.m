%% Load transformed data

clear
close all

load('Sacramento/transformedSacDataV2.mat')

%% cleaning up pass 527
% parts of the Sac are outside the swath of pass 526, so there are
% bad/missing values for many fields. Nodes outside the swath can be
% identified by having NaN value for nWidth, so all fields with nWidth =
% NaN are effectively deleted here. Set to NaN for array size continuity.

tic

for i = 1:length(truth)
    iNoWidth = isnan(truth(i).nWidth) | truth(i).geoHeight == -9999;

    simulated(i).reach(iNoWidth) = NaN;
    simulated(i).node(iNoWidth) = NaN;
    simulated(i).easting(iNoWidth) = NaN;
    simulated(i).northing(iNoWidth) = NaN;
    simulated(i).nHeight(iNoWidth) = NaN;
    simulated(i).nWidth(iNoWidth) = NaN;
    simulated(i).geoHeight(iNoWidth) = NaN;
    simulated(i).sCoord(iNoWidth) = NaN; 
    simulated(i).nCoord(iNoWidth) = NaN;


    truth(i).reach(iNoWidth) = NaN;
    truth(i).node(iNoWidth) = NaN;
    truth(i).easting(iNoWidth) = NaN;
    truth(i).northing(iNoWidth) = NaN;
    truth(i).nHeight(iNoWidth) = NaN;
    truth(i).nWidth(iNoWidth) = NaN;
    truth(i).geoHeight(iNoWidth) = NaN;
    truth(i).sCoord(iNoWidth) = NaN; 
    truth(i).nCoord(iNoWidth) = NaN;
end

nTotProf = length(simulated);
RMSEslmSUM = zeros(nTotProf,1);
MAEslmSUM = zeros(nTotProf,1);
RMSEsimAvgSUM = zeros(nTotProf,1);
nCombRec = zeros(nTotProf,1);

allSimulated = simulated;
allTruth = truth;


%% loop through number of observations
parfor numProfiles = 3:nTotProf
    disp(numProfiles)
    maxComb = 250;
    
    %generate all possible combinations of the profiles for current
    %numProfiles.
    allComb = combnk((1:nTotProf),numProfiles);
    
    if numProfiles == 1
        iValidComb = sum(mod(allComb,2),2)==1;
    else
    %only keep combinations that have >= 2 observations everywhere.
        iValidComb = sum((allComb<10),2)>=2;
        
    end
    allComb = allComb(iValidComb,:);
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


        % try some clever averaging

        % [simAvg, simIncl,S] = multiProfAvg(simulated);
        % [truthAvg, truthIncl] = multiProfAvg(truth,simIncl);
        if numProfiles ~= 1
            simAvg= nodeAvg2(simulated,0.5);
            truthAvg = nodeAvg2(truth,0.5);
%             simAvg = nodeAvg3(simulated,'geoHeight');
%             truthAvg = nodeAvg3(truth,'geoHeight');
        else
            simAvg = simulated;
            truthAvg = truth;
            wSim = ones(718,1);
            wTruth = ones(718,1);
        end

        % SLM toolbox 
        prescription = slmset('Decreasing','on');
        notNaN = ~isnan(simAvg.sCoord) & ~isnan(simAvg.geoHeight); %slm removes these
% 
%         if length(simulated)>2
%             prescription.Weights = simAvg.normWeight(notNaN);
%         end

        [~,iKnots] = slidePeaks(simAvg.nWidth(notNaN),0.10,0);

        locKnots = simAvg.sCoord(iKnots); %s-coordinate location for knots
        locKnots(isnan(locKnots)) = []; %delete NaNs
        locKnots = [min(simAvg.sCoord); locKnots; max(simAvg.sCoord)];

        prescription.Knots = locKnots;
        slm = slmengine(simAvg.sCoord(notNaN), simAvg.geoHeight(notNaN), prescription);
        slmProf = slmeval(slm.x,slm);  %evaluate
        knotZ = slmeval(slm.knots,slm);
        
%         
%         figure()
%         plot(simAvg.sCoord,simAvg.geoHeight,'o')
%         hold on
%         plot(slm.x,slmProf,'r','Linewidth',2)
%         if sqrt(nanmean((slmProf - truthAvg.geoHeight(notNaN)).^2)) < 2
%         end
        
        % errors
        RMSEslmSUM(numProfiles) = RMSEslmSUM(numProfiles) + sqrt(nanmean((slmProf - truthAvg.geoHeight(notNaN)).^2));
        MAEslmSUM(numProfiles,1) = MAEslmSUM(numProfiles,1) + nanmean(abs(slmProf - truthAvg.geoHeight(notNaN)));

        RMSEsimAvgSUM(numProfiles) =  RMSEsimAvgSUM(numProfiles) + sqrt(nanmean((simAvg.geoHeight - truthAvg.geoHeight).^2));

    end
    RMSEslm(numProfiles) = RMSEslmSUM(numProfiles) / nComb;
    MAEslm(numProfiles,1) = MAEslmSUM(numProfiles,1) /nComb;
    RMSEsimAvg(numProfiles) =  RMSEsimAvgSUM(numProfiles) /nComb;
end

toc
%% plots 
close all
plot(RMSEslm,'k','Linewidth',2)
hold on
ylabel('RMSE (m)')
xlabel('# of SWOT observations')
title('Average profile errors')
% hold off
% 



% 
% figure()
% histogram(simAvg.geoHeight(notNaN) - truthAvg.geoHeight)
% hold on
% histogram(slmProf - truthAvg.geoHeight)
% 


