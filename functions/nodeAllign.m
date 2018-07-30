function [profOut] = nodeAllign(prof, nodeField)
% Alligns all non-singeton fields from structure entries in prof into
% matricies and returns them in a structure with same field names. Optional
% input nodeField (default = 'node') specifies which field is used to
% allign the profiles.
%
% alligned = nodeAllign(simulated);
% alligned = nodeAllign(simulated,'fieldName');

if nargin < 2
    nodeField = 'node';
end

f = fields(prof);

%find required array dimensions
dim(2) = length(prof); %num profiles
minNode = Inf;
maxNode = -Inf;

%find min and max nodes
for i = 1:dim(2)
    minCheck = min(prof(i).(nodeField));
    maxCheck = max(prof(i).(nodeField));
    
    if  minCheck < minNode
        minNode = minCheck;
    end
    if maxCheck > maxNode
        maxNode = maxCheck;
    end
end

dim(1) = maxNode - minNode + 1; %range of nodes across all profiles

%initialize matricies
for i = 1:length(f)
    if size(prof(1).(f{i}),1) ~= 1
        profOut.(f{i})=nan(dim);
    end
end

fOut = fields(profOut);

for i = 1:dim(2)
    for j = 1:length(fOut)
        iNode = prof(i).node - minNode + 1;

        notNaN = ~isnan(iNode);
        iNode = iNode(notNaN);

        profOut.(fOut{j})(iNode,i) = prof(i).(fOut{j})(notNaN);
    end
end

profOut.nodeVec = (minNode:maxNode)';


