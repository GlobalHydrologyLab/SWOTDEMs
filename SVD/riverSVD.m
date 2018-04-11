%% SVD_multiSect_varReach.m Extending idea of SVD application to reduce
% errors in each profile individually to work for full river profiles
% instead of subsections. Reaches are defined by equiReach.m, which defines
% reaches only based on making them similar in length. 

%--------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------
% - look into removing erroneous data from Sacramento orbit 527 before SVD.
%       The large (and obvious) errors can dominate the decomposition
%       results and make comparisons of errors kind of disingenuous.
% x-- Hard-coded removal for now. I think this is fine.
%
%--------------------------------------------------------------------------

% clear
clearvars -except SVDStats SIMStats smoothStats opts
close all

opts.targetRL = 10;
opts.sectMin = 3;
opts.rmMean = 1;
opts.iSV = [1,3]; 
opts.maxDiff = 0.01; %set high to 'turn off' constraint.

% river = 'Sacramento';
% river = 'Po';
river = 'Tanana';

switch river
    case 'Sacramento'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Sacramento/SacDataV4.mat')
        zField = 'geoHeight';
        % % hard-coded removal of far range data from pass 527
        for i = 10:17
            simulated(i).geoHeight(339:487) = NaN;
        end
        simulated(6) = [];
        truth(6) = [];
        
    case 'Po'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Po/transformedPo_3Pass.mat')
        zField = 'nHeight';
        simulated(1:17) = trimFields(simulated(1:17),205:400);
        truth(1:17) = trimFields(truth(1:17),205:400);

    case 'Tanana'
        load('/Users/Ted/Documents/MATLAB/SWOTDEMs/Tanana/transformedTananaData.mat')
        zField = 'nHeight';
        simulated(9:12) = trimFields(simulated(9:12),250:675);
        truth(9:12) = trimFields(truth(9:12),250:675);
        simulated = trimFields(simulated,25:600);
        truth = trimFields(truth,25:600);
        opts.maxDiff = opts.maxDiff / 2; %100m node spacing
end

% clearvars -except simulated truth zField opts

%Group data into matrices without gaps, intelligently deleting data so that
%all sections have >= opts.sectMin rows
nProf = length(simulated);
simAllign = nodeAllign(simulated);
% zAll = simAllign.(zField);

%trim truth to sim extent
nodeRng = [min(simAllign.node(1,:)), max(simAllign.node(end,:))];
truth = trimFields(truth,nodeRng);

truthAllign = nodeAllign(truth);

%subSectByObs choses the rectangular matrices for svd
[section, zAll] = subsectByObs(simAllign.(zField),opts.sectMin);

%remove bias
observedBy = ~isnan(zAll);
truthAllign.(zField)(~observedBy) = NaN; %mask
zAll = zAll - nanmean(zAll - truthAllign.(zField),1);


%init. matrices for storing section data
dim = size(simAllign.sCoord);
z2All = nan(dim);
sAll = nan(dim);

for r = min(section):max(section)
    inSect = section == r;
    
    simReach = trimFields(simAllign,inSect);
    truthReach = trimFields(truthAllign,inSect);
    
    %assemble full rectangular matrices of s,z data 
    z = zAll(inSect,:);
%     z = truthAllign.(zField)(inSect,:); %testing svd on noise-free data
    [z, delCol] = nanRows(z,2);
    s = simAllign.sCoord(inSect,~delCol);

    %replace any missing s-values with mean of that node.
    [mMiss, nMiss] = find(isnan(s));
    for i = 1:length(mMiss)
        s(mMiss(i),nMiss(i)) = nanmean(s(mMiss(i),:));
    end
    
    if opts.rmMean
        %remove mean from each node.
        mz = nanmean(z,2);
        zresid = z - mz;
        [U,S,V] = svd(zresid);
    else
        [U,S,V] = svd(z);
    end
    

    %now modify S, removing smaller components.
%     z2 = SVRecomp(U,S,V,1:iSV(r));
    z2 = SVRecomp(U,S,V,opts.iSV);
    
    %recombine
    if opts.rmMean
        z2 = z2 + mz;
    end
  
    % join section data for later comparison  
    z2All(inSect,~delCol) = z2;
    sAll(inSect,~delCol) = s;
end

 
%constrain profiles
for i = 1:size(z2All,2)
%     weight(:,i) = movstd(zAll(:,i),7);
%     z2All(:,i) = slopeConstrain(z2All(:,i),opts.maxDiff, weight(:,i), 0);

    z2All(:,i) = slopeConstrain(z2All(:,i),opts.maxDiff);
end

% avgProf = slopeConstrain(z2All,opts.maxDiff);


skm = nanmean(sAll,2)/1000;

zErr = zAll - truthAllign.(zField);
z2Err = z2All - truthAllign.(zField);


simStats = reachStats(skm,zAll,truthAllign.(zField),opts.targetRL);
svdStats = reachStats(skm,z2All,truthAllign.(zField),opts.targetRL);
% avgStats = reachStats(skm,avgProf,nanmean(truthAllign.(zField),2), opts.targetRL);
%% 

%Renato's gaussian smoothing.
sigma = opts.targetRL/5;
for i = 1:size(zAll,2)
    [smoothProfs(:,i),~] = GaussianAveraging(skm,zAll(:,i), ... 
        simAllign.nWidth(:,i),opts.targetRL,sigma);
end
zSErr = smoothProfs - truthAllign.(zField);

smoothStats.(river) = reachStats(skm,smoothProfs,truthAllign.(zField),opts.targetRL);
smoothStats.(river).s = skm;
smoothStats.(river).z = smoothProfs;
smoothStats.(river).zErr = zSErr;

SVDStats.(river) = svdStats;
SVDStats.(river).s = skm;
SVDStats.(river).z = z2All;
SVDStats.(river).zErr = z2Err;
SVDStats.(river).x = simAllign.easting;
SVDStats.(river).y = simAllign.northing;
SVDStats.(river).meanProf = slopeConstrain(z2All,opts.maxDiff);

SIMStats.(river) = simStats;
SIMStats.(river).s = skm;
SIMStats.(river).z = zAll;
SIMStats.(river).zErr = zErr;
SIMStats.(river).x = simAllign.easting;
SIMStats.(river).y = simAllign.northing;

%% 
%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------

%rothko section plot
% figure()
% imagesc(~isnan(zAll) .* section)
% xlabel('Profile')
% ylabel('Node')
% c = lines;
% c = c(1:max(section),:);
% colormap([1 1 1; c])

% %coverage/elevation plot
% figure()
% imAlpha=ones(size(zAll));
% imAlpha(isnan(zAll))=0;
% imagesc(zAll,'AlphaData',imAlpha);
% colormap(brewermap(64,'YlGnBu'))
% set(gca,'color',0*[1 1 1]);

%singular values
figure()
bar(diag(S))
xlabel('Singular Value Number')
ylabel('Singular Value Magnitude')
% title('Singular Values of Elevation Residuals')
title(['Singular Values - ' sprintf(river)]);
% set(gca,'YScale','log');
set(gcf,'Position',[1000 987 862 351]);

%original and approx profiles
handle = figure();
subplot(2,1,1);
plot(skm,zAll)
% plot(skm,truthAllign.(zField))
hold on
% plot(skm,truthAllign.(zField)(:,3),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Original Simulation')

subplot(2,1,2);
plot(skm,z2All)
hold on
% plot(skm,truthAllign.(zField)(:,3),'k','Linewidth',2)
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Low-Rank Approximation')

% for r = reaches
%     ir = find(reachVec == r,1,'last');
%     xr(r) = skm(ir);
%     yr(r) = truthAllign.(zField)(ir,3);
% end
% subplot(2,1,1)
% scatter(xr,yr, 100, reaches,'filled')
% subplot(2,1,2)
% scatter(xr,yr, 100, reaches,'filled')

allAxes = findobj(handle, 'type', 'axes');
linkaxes(allAxes);
set(gcf,'Units','normalized','Position',[0.013672 0.013194 0.59922 0.91528])

% %epdf of all node errors
% figure()
% ksdensity(reshape(zErr,[],1))
% hold on
% ksdensity(reshape(z2Err,[],1))
% xlabel('Elevation Error (m)')
% title('EPDF of Node Errors')
% legend('Original Data','Low Rank')

figure()
subplot(2,2,1)
scatter(simStats.slopeMAE, svdStats.slopeMAE,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated MAE(cm/km)')
ylabel('LRA MAE(cm/km)')
title('MAE')
% 
subplot(2,2,2)
scatter(simStats.slopeRMSE, svdStats.slopeRMSE,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated RMSE(cm/km)')
ylabel('LRA RMSE(cm/km)')
title('RMSE')

subplot(2,2,3)
ksdensity(simStats.slopeErr(:))
hold on
ksdensity(svdStats.slopeErr(:))
legend('original','LRA')

svdStats.slopePctChange = (nanmean(svdStats.slopeMAE)-nanmean(simStats.slopeMAE))./nanmean(simStats.slopeMAE) * 100;
subplot(2,2,4)
text('Units','normalized','position',[0.05 0.8],'String', ... 
    ['Sim: ' num2str(nanmean(simStats.slopeMAE)) ' cm/km'], 'FontSize',24)
text('Units','normalized','position',[0.05 0.6],'String', ... 
    ['LRA: ' num2str(nanmean(svdStats.slopeMAE)) ' cm/km'], 'FontSize',24)
text('Units','normalized','position',[0.05 0.4],'String', ... 
    ['% change: ' num2str(svdStats.slopePctChange)], 'FontSize',24)
set(gca,'visible','off')
SVDStats.(river).slopePctChange = svdStats.slopePctChange;
% subplot(2,2,4)
% scatter(simStats.slopeRRMSE, svdStats.slopeRRMSE,[],nCol,'filled')
% scatter1to1(gca,'origin');
% xlabel('Simulated RRMSE(%)')
% ylabel('LRA RRMSE()')
% title('RRMSE')

% set(gcf,'Units','Normalized','Position',[0.58008 0.21389 0.37461 0.55833])

% reach elev. errors
% simStats.reachZMAE = nanmean(abs(simStats.reachZErr),1);
% svdStats.reachZMAE = nanmean(abs(svdStats.reachZErr),1);
% 
% simStats.reachZRMSE = sqrt(nanmean(simStats.reachZErr.^2,1));
% svdStats.reachZRMSE = sqrt(nanmean(svdStats.reachZErr.^2,1));
% 
% figure()
% plot(simStats.reachZRMSE', svdStats.reachZRMSE','k.','MarkerSize',20)
% % plot(abs(simStats.reachZErr), abs(svdStats.reachZErr),'ko')
% scatter1to1(gca);
% title('Reach-level elevation errors')
% xlabel('Simulated (m)')
% ylabel('LRA (m)')


% figure
% scatter(simStats.slopeRMSE, svdStats.slopeRMSE,'k','filled')
% scatter1to1(gca,'origin');
% xlabel('Simulated RMSE(cm/km)')
% ylabel('LRA RMSE(cm/km)')
% title(['Reach Slope Errors - ' river])
%--------------------------------------------------------------------------


%% gif
% 
% fileName = '/Users/Ted/Documents/MastersCommittee/firstMeeting/figs/method/16-1.gif';
% delete(fileName)
% figure()
% set(gcf,'Units','Normalized','Position',[0.319 0.495 0.352 0.25])
% plot(nanmean(s,2)/1000,SVRecomp(U,S,V,1:16)+mz)
% text('Units','normalized','position',[0.6 0.8],'String','rank: 16', 'FontSize',24)
% xlim = [16 26]; ylim = [28 35];
% set(gca,'XLim',xlim, 'YLim',ylim)
% box off
% xlabel('Flow Distance (km)')
% ylabel('Elevation (m)')
% gif(fileName,'DelayTime',0.75,'frame',gcf)
% 
% 
% for i = flip(1:15)
% plot(nanmean(s,2)/1000,SVRecomp(U,S,V,1:i)+mz)
% box off
% text('Units','normalized','position',[0.6 0.8],'String',['rank: ' num2str(i)], 'FontSize',24)
% set(gca,'XLim',xlim, 'YLim',ylim)
% xlabel('Flow Distance (km)')
% ylabel('Elevation (m)')
% gif
% end