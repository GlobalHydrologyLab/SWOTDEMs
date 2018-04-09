%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Renato.

clear

%find shapefiles in the directory
k = dir('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SWOTSimData/shapefiles/*.shp'); 
fileName = {k.name}';


isNodeFile = contains(fileName,'Node');
isTruthFile = contains(fileName,'Truth');


%all simulator data is processed to the same nodes, so instead of
%transforming each set of x,y coords, i simply replace with a previously
%transformed set in WGS84 UTM 10N.
load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SWOTSimData/eastNorthNodesUTM.mat')


% import shapfiles
sIndex=1;
tIndex=1;

for i = 1 : 1%length(fileName)
    
    shapefile = shaperead( [k(i).folder '/' fileName{i}] );
    
    %save node simulator files in struct 'truth'
    if isNodeFile(i) && ~isTruthFile(i) 
        
        simulated(sIndex).name = char(fileName(i));
        simulated(sIndex).reach = [shapefile.Reach_ID]';
        simulated(sIndex).node = [shapefile.Node_ID]';
        simulated(sIndex).easting = east;
        simulated(sIndex).northing = north;
        simulated(sIndex).nHeight = [shapefile.N_Hght]';
        simulated(sIndex).nWidth = [shapefile.N_width]';
        simulated(sIndex).geoHeight = [shapefile.Geo_hght]';
        simulated(sIndex).nObs = [shapefile.nobs]';
        
        sIndex = sIndex + 1;
    end
    
    
    %save node truth files in struct 'truth'
    if isNodeFile(i) && isTruthFile(i) 
        
        truth(tIndex).name = char(fileName(i));
        truth(tIndex).reach = [shapefile.Reach_ID]';
        truth(tIndex).node = [shapefile.Node_ID]';
        truth(tIndex).easting = east;
        truth(tIndex).northing = north;
        truth(tIndex).nHeight = [shapefile.N_Hght]';
        truth(tIndex).nWidth = [shapefile.N_width]';
        truth(tIndex).geoHeight = [shapefile.Geo_hght]';
        truth(tIndex).nObs = [shapefile.nobs]';
        
        tIndex = tIndex + 1;
    end
    
end

% clearvars -except simulated truth