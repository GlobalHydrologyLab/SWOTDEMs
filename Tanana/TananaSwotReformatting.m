%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Rui.

clear

%find shapefiles in the directory
k = dir('/Users/Ted/Documents/Tanana/SWOTSimulator/TananaRiverObs');
fileName = {k.name}';

isNodeFile = contains(fileName,'node');

% [N,refvec] = egm96geoid(1,[64.562 64.796],[-149.04 -147.94]);

sIdx=1;

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        csv = importdata( [k(i).folder '/' fileName{i}] );
            
        %three profiles have 0 elevations in a location...? This filters
        %them.
        hasElev = csv.data(:,13) ~= 0;
          
        simulated(sIdx).name = fileName{i}(1:end-4);
        simulated(sIdx).reach = flip(csv.data(hasElev,1));
        simulated(sIdx).node = 870 - flip(csv.data(hasElev,20)); %reverse the numerical order as well
        simulated(sIdx).easting = flip(csv.data(hasElev,22));
        simulated(sIdx).northing = flip(csv.data(hasElev,23));
        simulated(sIdx).nHeight = flip(csv.data(hasElev,13));
%         simulated(sIdx).geoHeightStd = flip(csv.data(hasElev,14));
        simulated(sIdx).nWidth = flip(csv.data(hasElev,8));
        simulated(sIdx).nWidthStd = flip(csv.data(hasElev,9));
        simulated(sIdx).nObs = flip(csv.data(hasElev,6));
        simulated(sIdx).lat = flip(csv.data(hasElev,2));
        simulated(sIdx).lon = flip(csv.data(hasElev,3));
       
        simulated(sIdx).geoHeight = simulated(sIdx).nHeight;

        sIdx = sIdx + 1;
   
    end
    
end


%% truth data
% sample input DEMs at the centerline locations

%reorder truth data so that we can directly compare with the simulator data
%by indices. 

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
    'elevs_69']; %from Elizabeth

k = dir('/Users/Ted/Documents/Tanana/SWOTSimulator/TananaRiverObs/*.tif');

for i = 1:length(k)
    p = [k(i).folder '/' k(i).name];
    tempName = k(i).name(1:end-4);
    
    s2 = cellstr(char(tempName));
    simID(i) = find(strcmp(simDataOrder,s2));
    
    cx = simulated(simID(i)).easting;
    cy = simulated(simID(i)).northing;
    
    prof = geotiffinterp(p, cx, cy,'nearest','buffer',11);
    
%     geoidCorrection = ltln2val(N,refvec, ... 
%         simulated(simID(i)).lat,simulated(simID(i)).lon,'bicubic');
%     prof = prof + geoidCorrection; %calculated ellipsoid correction
    %----------------------------------------------------------------------
%     prof = prof + 11.6; %ellipsoid correction per Elizabeth
    %----------------------------------------------------------------------
    
    
%     %masked areas (where centerline is outside river) have random large
%     %negative numbers, set to NaN here.
    
    t(simID(i)).name = tempName;
    t(simID(i)).simName = simulated(simID(i)).name;
    t(simID(i)).node = simulated(simID(i)).node;
    t(simID(i)).easting = cx;
    t(simID(i)).northing = cy;
    t(simID(i)).geoHeight = prof;
    
end

truth = t;

%%
simAllign = nodeAllign(simulated);
avgCenterline = [nanmean([simAllign.easting],2), ... 
    nanmean([simAllign.northing],2)];
transParam = [1 3 5 length(avgCenterline) 150]';

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

%% 

% clearvars -except simulated truth
% save('./Tanana/TananaSimData.mat')













