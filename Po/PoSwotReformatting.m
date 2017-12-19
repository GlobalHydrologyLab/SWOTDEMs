%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear


%find shapefiles in the directory
k = dir('/Users/Ted/Documents/Po/Pass532and461');
fileName = {k.name}';

isNodeFile = contains(fileName,'node');
isTruthFile = contains(fileName,'gdem');

% import shapfiles
sIndex=1;
tIndex=1;

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        
        shapefile = shaperead( [k(i).folder '/' fileName{i} '/' fileName{i} '.shp'] );
        
        %UTM transformation params
        if ~exist('utmstruct','var')
            poZone = utmzone(mean(shapefile(1).lat,'omitnan'), ... 
            mean(shapefile(1).lon,'omitnan'));
            
            utmstruct = defaultm('utm');
            utmstruct.zone = poZone;
            utmstruct.geoid = wgs84Ellipsoid;
            utmstruct = defaultm(utmstruct);
        end
        [east,north] = mfwdtran(utmstruct,[shapefile.Y]',[shapefile.X]');
        
        if ~isTruthFile(i) 
            simulated(sIndex).name = char(fileName(i));
            simulated(sIndex).reach = [shapefile.reach_indx]';
            simulated(sIndex).node = [str2num(char(shapefile.node_indx))];
            simulated(sIndex).easting = east;
            simulated(sIndex).northing = north;
            simulated(sIndex).nHeight = [shapefile.h_n_ave]';
            simulated(sIndex).nHeightStd = [shapefile.h_n_std]';
            simulated(sIndex).nWidth = [shapefile.w_ptp]';
            simulated(sIndex).nObs = [str2num(char(shapefile.nobs))];
            simulated(sIndex).lat = [shapefile.Y]';
            simulated(sIndex).lon = [shapefile.X]';
            
            sIndex = sIndex + 1;

        else
            truth(tIndex).name = char(fileName(i));
            truth(tIndex).reach = [shapefile.reach_indx]';
            truth(tIndex).node = [str2num(char(shapefile.node_indx))];
            truth(tIndex).easting = east;
            truth(tIndex).northing = north;
            truth(tIndex).nHeight = [shapefile.h_n_ave]';
            truth(tIndex).nHeightStd = [shapefile.h_n_std]';
            truth(tIndex).nWidth = [shapefile.w_ptp]';
            truth(tIndex).nObs = [str2num(char(shapefile.nobs))];
            truth(tIndex).lat = [shapefile.Y]';
            truth(tIndex).lon = [shapefile.X]';

            tIndex = tIndex + 1;
        end
    end
    
end


% % % % clearvars -except simulated truth
% % % % save('Po/PoSimulatorDataV2.mat')














