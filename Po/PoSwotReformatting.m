%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear


%find shapefiles in the directory
k = dir('/Users/Ted/Documents/Po/Ted_Po/*/');
fileName = {k.name}';

isNodeFile = contains(fileName,'node');
isTruthFile = contains(fileName,'gdem');

% import shapfiles
sIdx=1;
tIdx=1;

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        
        shapefile = shaperead( [k(i).folder '/' fileName{i} '/' fileName{i} '.shp'] );
        
        %UTM transformation params
        if ~exist('utmstruct','var')
            poZone = utmzone(nanmean(shapefile(1).lat), ... 
            nanmean(shapefile(1).lon));
            
            utmstruct = defaultm('utm');
            utmstruct.zone = poZone;
            utmstruct.geoid = wgs84Ellipsoid;
            utmstruct = defaultm(utmstruct);
        end
        [east,north] = mfwdtran(utmstruct,[shapefile.Y]',[shapefile.X]');
        
        if ~isTruthFile(i) 
            simulated(sIdx).name = char(fileName(i));
            simulated(sIdx).reach = [shapefile.reach_indx]';
            simulated(sIdx).node = [shapefile.node_indx]';
            simulated(sIdx).easting = east;
            simulated(sIdx).northing = north;
            simulated(sIdx).nHeight = [shapefile.h_n_ave]';
            simulated(sIdx).nHeightStd = [shapefile.h_n_std]';
            simulated(sIdx).nWidth = [shapefile.w_ptp]';
            simulated(sIdx).nObs = [shapefile.nobs]';
            simulated(sIdx).lat = [shapefile.Y]';
            simulated(sIdx).lon = [shapefile.X]';
            
            sIdx = sIdx + 1;

        else
            truth(tIdx).name = char(fileName(i));
            truth(tIdx).reach = [shapefile.reach_indx]';
            
            %weird workaround for some files being interpreted with
            %character values for nodes and nObs. Doesn't happen with
            %simulated files.
            if ischar([shapefile.node_indx])
                for j = 1:length(shapefile)
                    truth(tIdx).node(j,1) = str2double(shapefile(j).node_indx);
                end
            else
                truth(tIdx).node = [shapefile.node_indx]';
            end
            
            truth(tIdx).easting = east;
            truth(tIdx).northing = north;
            truth(tIdx).nHeight = [shapefile.h_n_ave]';
            truth(tIdx).nHeightStd = [shapefile.h_n_std]';
            truth(tIdx).nWidth = [shapefile.w_ptp]';

            if ischar([shapefile.nobs])
                for j = 1:length(shapefile)
                    truth(tIdx).nObs(j,1) = str2double(shapefile(j).nobs);
                end
            else
                truth(tIdx).nObs = [shapefile.nobs]';
            end
            truth(tIdx).lat = [shapefile.Y]';
            truth(tIdx).lon = [shapefile.X]';

            tIdx = tIdx + 1;
        end
    end
    
end
%%

% now we have a problem where truth files and simulator files do not match
% up. Left and right simulator outputs in pass 461 are separate files,
% whereas truth inputs are one profile for that overpass. This is the best
% way I can think to pair and pare the files.

for i = 1:numel(truth)
    iPD = regexp(truth(i).name,'[1234567890]'); %find any numbers in name
    tPassDay(i,:) = str2double(truth(i).name(iPD));
end

%now find matching truth file for each simulator file and trim to extent.
for i = 1:numel(simulated)
    iPD = regexp(simulated(i).name,'[1234567890]');
    sPassDay(i,:) = str2double(simulated(i).name(iPD));
    iT = find(sPassDay(i) == tPassDay);
    nodes = simulated(i).node;
    
    truthLR(i) = trimFields(truth(iT),nodes);
end

truth = truthLR;

% % % % % % clearvars -except simulated truth
% % % % % % save('Po/PoSim_3Pass.mat')














