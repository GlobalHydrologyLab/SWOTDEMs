function [profOut, observedBy] = nodeAvg3_1(prof,zField,profMask)
%Node average3, same purpose as previous iterations, but different
%conceptual approach. Want to test for signficant difference between input
%profiles over sections. Sections defined where all profiles either exist
%or don't exist everywhere, i.e. no partial profiles within sections. Must
%be able to handle profiles of different lengths and with different (but
%similar) node locations, provided there are consistent nodeID's for all
%profiles. 

%profMask is a logical array of profiles to include, intended to limit
%truth profiles to those that were available from the simulated data.

%11/6/17 solves same as nodeAvg3 but more generally for all
%non-singleton fields. nodeAvg3 had hardcoded fields.

%   TO DO:
% X subsectByObs.m - identify suitable sections as defined above. 
% X remove mean z of each profile (need to think about this more)
%   stackProfiles.m -  written but not implemented.
% - test (paired rank?) profile sections for similarity
% - average profiles that are sufficiently similar
% ? add mean z of each profile back.
% - join sections

if length(prof)==1
    warning('Only one profile given, returning input')
    profOut = prof;
    profOut.normWeight = ones(length(prof.node),1);
    observedBy = ~isnan(prof.(zField));
    return
end

profAllign = nodeAllign(prof);
f = fields(profAllign);
numFields = length(f);
dim = size(profAllign.node);

%if no mask supplied, create all true mask.
if nargin < 3
    profMask = ones(dim);
end

%apply mask
for i = 1:numFields
    profAllign.(f{i})(~profMask) = NaN;
end

z = profAllign.(zField);
avgZ = nanmean(z,2); %simple for testing

%find sections of common observations
observedBy = ~isnan(z);
 
S = nanstd(z,0,2);
w = 1./S;
 
%weight assignment breaks when S = 0 -> w = Inf, which only happens when
%there is just 1 observation at that node. assign weight 1 in that case.
w(w==Inf) = 1;

%normalize such that mean(wNorm) = 1
wNorm = length(z) .* w ./ sum(w);

%assemble return profile
profOut.(zField) = avgZ;
% profOut.nWidth = nanmedian(profAllign.nWidth,2);
profOut.normWeight = wNorm;
profOut.node = profAllign.nodeVec;

fExist = fields(profOut);
for i = 1:numFields
    if ~any(strcmp(f{i},fExist))
        profOut.(f{i}) = nanmean(profAllign.(f{i}),2);
    end
end

end

