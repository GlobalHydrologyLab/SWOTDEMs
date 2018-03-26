function [out] = stackProfiles(prof, zField)
% removes mean elevation differences from profiles. Finds the longest
% profile for reference and removes difference in mean elevations between
% all other profiles and the reference profile over the section they
% overlap. note: ref must overlap in s coordinate with all other profiles,
% more overlap will reduce effect of noise on offset. Node IDs must be
% consistent between profiles, as allignment depends on nodes.

d = nodeAllign(prof);


hasObs = ~isnan(d.node);
[~,longestID] = max(sum(hasObs,1));

z = d.(zField);
z(z==0) = NaN;
s = d.sCoord;

refZ = z(:,longestID);% - nanmean(z(:,longestID));

for i = 1:length(prof)
    %only has result in overlapping sections
    meanDiff(i) = nanmean(z(:,i) - refZ);
    
    z0(:,i) = z(:,i) - meanDiff(i);
end

z0 = z0 + mean(meanDiff);

sV = [];
zV = [];
for i = 1:length(prof)
    hasData = ~isnan(s(:,i)) & ~isnan(z(:,i));
    sV = [sV; s(hasData,i)];
    zV = [zV; z0(hasData,i)];
end

[sV, sortOrder] = sort(sV);
zV = zV(sortOrder);

out.s = s;
out.z = z0;
out.sV = sV;
out.zV = zV;
out.nodeIDs = d.nodeVec;

end

