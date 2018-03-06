% minRL.m
%
% Ideas on solving the minimum reach length from LeFavour and Alsdorf more
% generally. At finer resolutions than their data, minimum slopes can be
% very small, leading to excessively large reach lengths.

clear
close all

load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/transformedSacDataV2.mat')
zField = 'geoHeight';


clearvars -except simulated truth zField

%resample data to larger node spacings and show change in RL

dx = 10:10:5000;

for p = 1:length(simulated)
    
    xMin = min(simulated(p).sCoord);
    xMax = max(simulated(p).sCoord);
    simulated(p).(zField)(simulated(p).(zField) == -9999) = NaN;
    truth(p).(zField)(truth(p).(zField) == -9999) = NaN;
    
    for i = 1:length(dx)
       d = dx(i);
       sRng = xMin:d:xMax;
       
       zHat = interp1(simulated(p).sCoord, simulated(p).(zField), ...
           sRng,'linear');
        
       tzHat = interp1(truth(p).sCoord, truth(p).(zField), ...
           sRng, 'linear');
       
       err = zHat - tzHat;
       sd = std(err,'omitnan');
       
       tSlope = min(abs(diff(tzHat) ./ diff(sRng)));
  
       
       RL(i,p) = 2*sd/tSlope;
       sMin(i,p) = tSlope;
       
       
    end


end

%%

figure()
plot(dx,RL./1000)
set(gca, 'YScale', 'log')

figure()
plot(dx,sMin)

