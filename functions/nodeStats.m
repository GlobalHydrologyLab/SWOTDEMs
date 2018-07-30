function [testStats] = nodeStats(skm,testZ,truthZ)
%%NODESTATS
%Takes in distance vector (skm), test (testZ) and truth (truthZ) elevation
%profiles and calculates errors, MAE, and RMSE. MAE and RMSE are calculated
%for each column (overpass) as well as total from all height measurements.
%
%testStats = NODESTATS(skm,testZ,truthZ);

if size(testZ,2) ~= size(truthZ,2)
    error(['Arguments testZ and truthZ must have the same number of' ... 
        ' profiles for comparison']);
end

testStats.s = skm;
testStats.z = testZ;
testStats.zErr = testZ - truthZ;
testStats.dailyMAE = nanmean(abs(testStats.zErr),1);
testStats.dailyRMSE = sqrt(nanmean((testStats.zErr).^2,1));
testStats.totMAE = nanmean(abs(testStats.zErr(:)));
testStats.totRMSE = sqrt(nanmean((testStats.zErr(:)).^2));
end

