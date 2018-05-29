%% SacSwotReformatting.m 
%Reads in all shapefiles in a directory and saves them as structures. Set
%up according to format of shapefiles of the Sac given by Renato.

clear

%find shapefiles in the directory
k = dir('./Sacramento/SWOTSimData/shapefiles/*.shp'); 
fileName = {k.name}';

isNodeFile = contains(fileName,'Node');
isTruthFile = contains(fileName,'Truth');

utmStruct = defaultm('utm'); 
utmStruct.zone = '10N';  
utmStruct.geoid = wgs84Ellipsoid;
utmStruct = defaultm(utmStruct);

% import shapfiles
sIndex=1;
tIndex=1;

load geoid

for i = 1 : length(fileName)
    
    if isNodeFile(i)
        shapefile = shaperead( [k(i).folder '/' fileName{i}] );
    end
    
    %save node simulator files in struct 'truth'
    if isNodeFile(i) && ~isTruthFile(i) 
        
        simulated(sIndex).name = char(fileName(i));
        simulated(sIndex).reach = [shapefile.Reach_ID]';
        simulated(sIndex).node = [shapefile.Node_ID]';
        simulated(sIndex).lat = [shapefile.Y]';
        simulated(sIndex).lon = [shapefile.X]';
        [simulated(sIndex).easting, simulated(sIndex).northing] = ...
            mfwdtran(utmStruct,simulated(sIndex).lat,simulated(sIndex).lon);
        simulated(sIndex).nHeight = [shapefile.N_Hght]';
        simulated(sIndex).nWidth = [shapefile.N_width]';
        simulated(sIndex).nObs = [shapefile.nobs]';
        
        missingData = simulated(sIndex).nHeight == -9999;
        simulated(sIndex).nHeight(missingData) = NaN;
        simulated(sIndex).nWidth(missingData) = NaN;
        
        simulated(sIndex).geoHeight = simulated(sIndex).nHeight - ... 
            mean(ltln2val(geoid, geoidrefvec, ...
            simulated(sIndex).lat, simulated(sIndex).lon));
        
        sIndex = sIndex + 1;
    end
    
    
    %save node truth files in struct 'truth'
    if isNodeFile(i) && isTruthFile(i) 
        
        truth(tIndex).name = char(fileName(i));
        truth(tIndex).reach = [shapefile.Reach_ID]';
        truth(tIndex).node = [shapefile.Node_ID]';
        truth(tIndex).lat = [shapefile.Y]';
        truth(tIndex).lon = [shapefile.X]';
        [truth(tIndex).easting, truth(tIndex).northing] = ...
            mfwdtran(utmStruct,truth(tIndex).lat,truth(tIndex).lon);
        truth(tIndex).nHeight = [shapefile.N_Hght]';
        truth(tIndex).nWidth = [shapefile.N_width]';
        truth(tIndex).nObs = [shapefile.nobs]';
        
        missingData = truth(tIndex).nHeight == -9999;
        truth(tIndex).nHeight(missingData) = NaN;
        truth(tIndex).nWidth(missingData) = NaN;
        
        truth(tIndex).geoHeight = truth(tIndex).nHeight - ... 
            mean(ltln2val(geoid, geoidrefvec, ...
            truth(tIndex).lat, truth(tIndex).lon));
        
        tIndex = tIndex + 1;
    end
    
end

trimRng = 1:715;
simulated = trimFields(simulated,trimRng);
truth = trimFields(truth,trimRng);

%calc s n coords
avgCenterline = [nanmean([simulated.easting],2), ... 
    nanmean([simulated.northing],2)];
transParam = [0 3 5 length(avgCenterline) 100]';

for i = 1:length(simulated)
   sn = xy2sn(avgCenterline, ... 
       [simulated(i).easting,simulated(i).northing], transParam);
   simulated(i).sCoord = sn(:,1);
   simulated(i).nCoord = sn(:,2);
   xy2snAscendCheck(sn(:,1));
   
   sn = xy2sn(avgCenterline, ... 
       [truth(i).easting, truth(i).northing], transParam);
   truth(i).sCoord = sn(:,1);
   truth(i).nCoord = sn(:,2);
   xy2snAscendCheck(sn(:,1));
end


clearvars -except simulated truth
save('./Sacramento/SacSimData.mat')



