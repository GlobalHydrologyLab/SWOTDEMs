%% TananaTruthData.m
%Samples input geotiffs at the node locations to get truth profiles.
%Probably should redo this riverObs.
clear
close all


load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/avgTananaCenterline.mat')
cx = centerline(:,1);
cy = centerline(:,2);

k = dir('/Users/Ted/Documents/Tanana/SWOTSimulator/TananaRiverObs/*.tif');

for i = 1:length(k)
    p = [k(i).folder '/' k(i).name];
    
    prof = geotiffinterp(p, cx, cy,'nearest','buffer',3,'show');
    
    %----------------------------------------------------------------------
    prof = prof + 11.6; %ellipsoid correction per Elizabeth
    %----------------------------------------------------------------------
    
    
%     %masked areas (where centerline is outside river) have random large
%     %negative numbers, set to NaN here.
    
    t(i).name = k(i).name(1:end-4);
    t(i).node = clNodes;
    t(i).easting = cx;
    t(i).northing = cy;
    t(i).nHeight = prof;
    
end

truth = t;

%% 
% reorder truth data so that we can directly compare with the simulator
% data by indecies.

% simDataOrder = ['elevs_2 ';
%     'elevs_23';
%     'elevs_44';
%     'elevs_65';
%     'elevs_5 ';
%     'elevs_26';
%     'elevs_47';
%     'elevs_68';
%     'elevs_6 ';
%     'elevs_27';
%     'elevs_48';
%     'elevs_69'];
% 
% 
% s2 = cellstr(char(t.name));
% 
% for i = 1:length(simDataOrder)
%     match(i) = strcmp(simDataOrder(i,:),s2);
%     
% end



% % % % % clearvars -except truth
% % % % % save('Tanana/InterpTruth.mat')