function [profiles] = trimFields(profiles, trimBy)
%% trimFields.m
% Trims all fields in all profiles containing more than one value according
% to nodeMask. Assumed that non-single value fields have length equal to
% number of nodes. trimBy is vector with length equal to number of nodes,
% where 0 marks deletion of value. Alternatively, trimBy can be the range
% of node IDs- in this case a binary mask will be created where
% profiles.nodes is within the range of trimBy. trimBy is not to be
% confused with the greasy hat, trilby.
%
% Function written to trim down simulator input data set to a desired
% extent. 11/1/17.

numProfiles = size(profiles,2);

for i = 1:numProfiles    
    fieldNames = fields(profiles(i));
    numFields = length(fieldNames);
    
    %if the trimBy is node range, create logical masks
    if ~islogical(trimBy)
        nodeMask = profiles(i).node >= min(trimBy) & ... 
            profiles(i).node <= max(trimBy);
    else
        nodeMask = trimBy;
    end
        
    for j = 1:numFields
       
        %make sure it's not a single value field (e.g. 'name').
        if size(profiles(i).(fieldNames{j}),1) ~= 1
            %what in(dex) tarnation
            profiles(i).(fieldNames{j})(~nodeMask,:) = [];    
        end
    end
end

end

