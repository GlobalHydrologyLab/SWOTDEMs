function [testStats] = reachStats(skm,testZ,truthZ,targetRL)
%%reachStats.m
%Takes in distance vector (skm), test (testZ) and truth (truthZ) elevation
%profiles and calculates slope differences with a moving window of length
%~targetRL. Does not resample points, so the window is likely not exactly
%targetRL long. skm and targetRL are the same units.
%
%testZ and truthZ profiles are arranged such that rows are locations along
%river profile, and columns are different days. skm, testZ, and truthZ,
%must all be the same number of rows, with each row representing the same
%location between variables. Output arg is arranged similarly, but less
%rows due to the window not sliding past the ends.
%
%testStats = reachStats(skm,testZ,truthZ,targetRL);

if size(skm,1) ~= size(testZ,1) || size(skm,1) ~= size(truthZ,1)
    error(['Arugments skm, testZ, and truthZ must all have the same '... 
        'number of rows.']);
elseif size(testZ,2) ~= size(truthZ,2)
    error(['Arguments testZ and truthZ must have the same number of' ... 
        ' profiles for comparison']);
end

[nNodes,nProf] = size(testZ);

meanNodeSpc = abs(nanmean(diff(skm)));
nodesPerReach = round(targetRL/meanNodeSpc);

nReaches = nNodes - nodesPerReach;

testStats.slope = nan(nReaches,nProf);
testStats.slopeErr = nan(nReaches,nProf);
testStats.reachAvgZErr = nan(nReaches,nProf);

for r = 1:nReaches
    
    reachIds = r : (r+nodesPerReach-1);
    
    for p = 1:nProf
        testZRP = testZ(reachIds,p);
        truthZRP = truthZ(reachIds,p);
        skmRP = skm(reachIds);
        
        [testZRP,nanIds] = nanRows(testZRP);
        truthZRP(nanIds) = [];
        skmRP(nanIds) = [];
        
        if sum(~isnan(testZRP)) >= 0.9*nodesPerReach
            fitTest = polyfit(skmRP,testZRP,1);
            fitTruth = polyfit(skmRP,truthZRP,1);
            
            testStats.slope(r,p) = fitTest(1)*100;
            testStats.slopeErr(r,p) = fitTest(1)*100 - fitTruth(1)*100;
            testStats.zErr(r,p) = mean(testZRP) - mean(truthZRP);
            
        end
        
    end
    
end

testStats.slopeMAE = nanmean(abs(testStats.slopeErr),1);
testStats.slopeRMSE = sqrt(nanmean(testStats.slopeErr.^2,1));
testStats.nodeMAE = nanmean(abs(testZ - truthZ),1);
testStats.nodeRMSE = sqrt(nanmean((testZ - truthZ).^2,1));

end

