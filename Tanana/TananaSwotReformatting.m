%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear

zField = 'nHeight';

%find shapefiles in the directory
k = dir('/Users/Ted/Documents/Tanana/SWOTSimulator/TananaRiverObs');
fileName = {k.name}';

isNodeFile = contains(fileName,'node');

sIndex=1;

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        csv = importdata( [k(i).folder '/' fileName{i}] );
            
        %three profiles have 0 elevations in a location...? This filters
        %them.
        hasElev = csv.data(:,13) ~= 0;
          
        simulated(sIndex).name = fileName{i}(1:end-4);
        simulated(sIndex).reach = flip(csv.data(hasElev,1));
        simulated(sIndex).node = flip(csv.data(hasElev,20));
        simulated(sIndex).easting = flip(csv.data(hasElev,22));
        simulated(sIndex).northing = flip(csv.data(hasElev,23));
        simulated(sIndex).nHeight = flip(csv.data(hasElev,13));
        simulated(sIndex).nHeightStd = flip(csv.data(hasElev,14));
        simulated(sIndex).nWidth = flip(csv.data(hasElev,8));
        simulated(sIndex).nObs = flip(csv.data(hasElev,6));
        simulated(sIndex).lat = flip(csv.data(hasElev,2));
        simulated(sIndex).lon = flip(csv.data(hasElev,3));

        sIndex = sIndex + 1;

        
    end
    
end

%% truth data
% sampled input DEMs at the centerline locations

%reorder truth data so that we can directly compare with the simulator data
%by indeces.

simDataOrder = ['elevs_2 ';
    'elevs_23';
    'elevs_44';
    'elevs_65';
    'elevs_5 ';
    'elevs_26';
    'elevs_47';
    'elevs_68';
    'elevs_6 ';
    'elevs_27';
    'elevs_48';
    'elevs_69'];

k = dir('/Users/Ted/Documents/Tanana/SWOTSimulator/TananaRiverObs/*.tif');

for i = 1:length(k)
    p = [k(i).folder '/' k(i).name];
    t(i).name = k(i).name(1:end-4);
    
    s2 = cellstr(char(t(i).name));
    simID(i) = find(strcmp(simDataOrder,s2));
    
    cx = simulated(simID(i)).easting;
    cy = simulated(simID(i)).northing;
    
    prof = geotiffinterp(p, cx, cy,'nearest','buffer',3,'show');
    
    %----------------------------------------------------------------------
    prof = prof + 11.6; %ellipsoid correction per Elizabeth
    %----------------------------------------------------------------------
    
    
%     %masked areas (where centerline is outside river) have random large
%     %negative numbers, set to NaN here.
    
    t(i).simName = simulated(simID(i)).name;
    t(i).node = simulated(simID(i)).node;
    t(i).easting = cx;
    t(i).northing = cy;
    t(i).nHeight = prof;
    
end

truth = t;

%% 

clearvars -except simulated truth zField
save('Tanana/TananaSimTruth.mat')














