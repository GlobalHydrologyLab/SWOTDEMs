%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear


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

zField = 'nHeight';

clearvars -except simulated truth zField
save('Tanana/TananaSimulatorData.mat')














