%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear


%find shapefiles in the directory
k = dir('/Users/Ted/Documents/Po/Ted_Po/*/');
fileName = {k.name}';

isNodeFile = contains(fileName,'node');
isTruthFile = contains(fileName,'gdem');

utmStruct = defaultm('utm'); 
utmStruct.zone = '32T';  
utmStruct.geoid = wgs84Ellipsoid;
utmStruct = defaultm(utmStruct);

% import shapfiles
sIdx=1;
tIdx=1;

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        
        shapefile = shaperead( [k(i).folder '/' fileName{i} '/' fileName{i} '.shp'] );
        
        
        if ~isTruthFile(i) 
            simulated(sIdx).name = char(fileName(i));
            simulated(sIdx).reach = [shapefile.reach_indx]';
            simulated(sIdx).node = [shapefile.node_indx]';
            simulated(sIdx).lat = [shapefile.Y]';
            simulated(sIdx).lon = [shapefile.X]';
            [simulated(sIdx).easting, simulated(sIdx).northing] = ...
                mfwdtran(utmStruct,simulated(sIdx).lat,simulated(sIdx).lon);
            simulated(sIdx).nHeight = [shapefile.h_n_ave]';
            simulated(sIdx).nHeightStd = [shapefile.h_n_std]';
            simulated(sIdx).nWidth = [shapefile.w_ptp]';
            simulated(sIdx).nObs = [shapefile.nobs]';

            
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
            
            truth(tIdx).lat = [shapefile.Y]';
            truth(tIdx).lon = [shapefile.X]';
            [truth(tIdx).easting, truth(tIdx).northing] = ...
                mfwdtran(utmStruct,truth(tIdx).lat,truth(tIdx).lon);
            truth(tIdx).nHeight = [shapefile.h_n_ave]';
            truth(tIdx).nHeightStd = [shapefile.h_n_std]';
            truth(tIdx).nWidth = [shapefile.w_ptp]';

            if ischar([shapefile.nobs]) %same weird character issues
                for j = 1:length(shapefile)
                    truth(tIdx).nObs(j,1) = str2double(shapefile(j).nobs);
                end
            else
                truth(tIdx).nObs = [shapefile.nobs]';
            end

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

%%

simAllign = nodeAllign(simulated);
avgCenterline = [nanmean([simAllign.easting],2), ... 
    nanmean([simAllign.northing],2)];
transParam = [1 3 5 length(avgCenterline) 200]';

for i = 1:length(simulated)
   [sn,~,clOut] = xy2sn(avgCenterline, ... 
       [simulated(i).easting,simulated(i).northing], transParam);
   simulated(i).sCoord = sn(:,1);
   simulated(i).nCoord = sn(:,2);
   
   sn = xy2sn(avgCenterline, ... 
       [truth(i).easting, truth(i).northing], transParam);
   truth(i).sCoord = sn(:,1);
   truth(i).nCoord = sn(:,2);
end


clearvars -except simulated truth
save('Po/PoSimData.mat')














