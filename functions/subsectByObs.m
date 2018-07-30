function [section, zArray] = subsectByObs(zArray, sectMin)
%% subsectByObs.m
% written for SWOT DEM project.
%
% Input is array with columns of river profiles, with rows according to
% river nodes. Unobserved node/profile combinations are identified as NaN.
% Function nodeAllign will create zArray from SWOTDEM structure.
%
% returns vector (length = numNodes) of section IDs.
%
% Function finds sections of the river that are observed by the same
% combination of profiles.
%
% optional merging of sections with less than sectMin nodes. Finds shortest
% section that is less than sectMin and merges with neighboring section
% which requires fewer deletions of observations. This is repeated until
% there are no remaining sections that are shorter than sectMin.
%
% Currently the merge option forces mergers until all sections are longer
% than sectMin, which can have undesireable results if sectMin is set high.

% Ted Langhorst

%%

dim = size(zArray);
observedBy = ~isnan(zArray);

% initial section classification. Checks pattern of observations against
% previous to determine if current node is new section.
pattern = observedBy(1,:);

section = zeros(dim(1),1);
section(1) = 1;

sectNum = 1;
for i = 2:dim(1)
    %check if node is observed by different profiles
    if any(observedBy(i,:) ~= pattern)
        sectNum = sectNum + 1;
        pattern = observedBy(i,:);
    end
    
    section(i) = sectNum;
end
 

% merge sections with few nodes, remove observations which make them
% different. 
if exist('sectMin','var')
    
    if sectMin > dim(1)
        error('sectMin must be <= number of rows in data')
    end
    
    nObs = sum((section==1:max(section))); %number of nodes in each section ID
    shortList = find(nObs < sectMin); %section IDs of short sections
    [~,shortSort] = sort(nObs(shortList));
    shortList = shortList(shortSort); %ascending length sections
    
    try
        toMerge = shortList(1);
    catch 
        toMerge = [];
    end
    
    %merge sections until there aren't any that are too small.
    while ~isempty(toMerge)
        nDelUp = Inf; %starting values
        nDelDown = Inf;
        
        %check number of deletions required for merger up
        if toMerge ~= min(section)
            bothSectUp = (section == toMerge | section == toMerge-1);
            obs = observedBy(bothSectUp, :);
            
            %column(s) with partial obs over potential merged section.
            badColUp = (sum(obs,1) ~= size(obs,1) & sum(obs,1) ~= 0);
            
            nDelUp = sum(sum(obs(:,badColUp)));
        end
        
        %check number of deletions required for merger down
        if toMerge ~= max(section)
            bothSectDown = (section == toMerge | section == toMerge+1);
            obs = observedBy(bothSectDown, :);
            
            %column(s) with partial obs over potential merged section.
            badColDown = (sum(obs,1) ~= size(obs,1) & sum(obs,1) ~= 0);
            
            nDelDown = sum(sum(obs(:,badColDown)));
        end
        
        %check self deletion case
        selfSect = section == toMerge;
        obs = observedBy(selfSect,:);
        nDelSelf = sum(sum(obs));
        
        
        %Do merger which requires least deletions. Delete necessary obs
        %from zArray, update observedBy, update section values to reflect
        %merger.
        mergeUp = nDelUp <= nDelDown & nDelUp < nDelSelf;
        mergeDown = nDelDown < nDelUp & nDelDown < nDelSelf;
        
        if mergeUp
            zArray(bothSectUp,badColUp) = NaN;
            observedBy(bothSectUp,badColUp) = 0;
            section(section>=toMerge) = section(section>=toMerge) - 1;
%             deleted(bothSectUp,badColUp) = 1;
        elseif mergeDown
            zArray(bothSectDown,badColDown) = NaN;
            observedBy(bothSectDown,badColDown) = 0;
            section(section>toMerge) = section(section>toMerge) - 1;
%             deleted(bothSectDown,badColDown) = 1;
        else
            zArray(selfSect,:) = NaN;
            observedBy(selfSect,:) = 0;
            section(section==toMerge) = NaN;
            section(section>toMerge) = section(section>toMerge) - 1;
%             deleted(selfSect,:) = 1;
        end
        
        %update list of short sections
        nObs = sum((section==1:max(section)),1); %number of nodes in each section ID
        shortList = find(nObs < sectMin); %section IDs of short sections
        [~,shortSort] = sort(nObs(shortList));
        shortList = shortList(shortSort); %ascending length sections

        
        try
            toMerge = shortList(1);
        catch 
            toMerge = [];
        end
    end
    
    %final check to see if there are any sections that can be merged. This
    %can happen when both sections are larger than the minimum section
    %length, as the merging loop will have stopped.
    i=1;
    while i < max(section)
        
        thisSect = observedBy(section==i,:);
        nextSect = observedBy(section==i+1,:);
        if all(thisSect(1,:) == nextSect(1,:))
            section(section>i) = section(section>i) - 1; 
        else
            %only iterate if a section is not removed. Otherwise a section
            %is skipped due to the reduced section number.
            i = i + 1;
        end
    end
    
end

end

