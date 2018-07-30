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
clearvars -except SVDStats SIMStats smoothStats rOpts
close all

opts.sectMin = 2;
opts.maxDiff = 0; %set high to 'turn off' constraint.
zField = 'geoHeight';

river = 'Sacramento';
% river = 'Po'; 
% river = 'Tanana';

switch river
    case 'Sacramento'
        load('./Sacramento/SacSimData.mat')
        % % hard-coded removal of far range data from pass 527
%         for i = 10:17
%             simulated(i).(zField)(339:487) = NaN;
%         end
        simulated(6) = []; %remove high discharge day
        truth(6) = [];
        groups = [ones(1,8) ones(1,8)+1]';
        
    case 'Po'
        load('./Po/PoSimData.mat')
        simulated(1:17) = trimFields(simulated(1:17),205:400);
        truth(1:17) = trimFields(truth(1:17),205:400);
        groups = [ones(1,17) ones(1,34)+1 ones(1,18)+2]';

    case 'Tanana'
        load('./Tanana/TananaSimData.mat') 
        simulated(9:12) = trimFields(simulated(9:12),250:675);
        truth(9:12) = trimFields(truth(9:12),250:675);
        simulated = trimFields(simulated,25:600);
        truth = trimFields(truth,25:600);
        opts.maxDiff = opts.maxDiff / 2; %100m node spacing
        groups = [ones(1,4) ones(1,4)+1 ones(1,4)+2]';
end

% clearvars -except simulated truth zField opts

%Group data into matrices without gaps, intelligently deleting data so that
%all sections have >= opts.sectMin rows
simAllign = nodeAllign(simulated);

%trim truth to sim extent
nodeRng = [min(simAllign.node(1,:)), max(simAllign.node(end,:))];
truth = trimFields(truth,nodeRng);
truthAllign = nodeAllign(truth);

%subSectByObs chooses the rectangular matrices for svd
[section, zAll] = subsectByObs(simAllign.(zField),opts.sectMin);
% [section, zAll] = subsectByObs(truthAllign.(zField),opts.sectMin);

%remove bias
observedBy = ~isnan(zAll);
truthAllign.(zField)(~observedBy) = NaN; %mask
zAll = zAll - nanmean(zAll - truthAllign.(zField),1);

%init. matrices for storing section data
dim = size(simAllign.sCoord);
z2All = nan(dim);
sAll = nan(dim);
opts.avgSV = 0;

for r = min(section):max(section)
    inSect = section == r;
    
    simReach = trimFields(simAllign,inSect);
    truthReach = trimFields(truthAllign,inSect);
    
    %assemble full rectangular matrices of s,z data 
    z = zAll(inSect,:);
    [z, delCol] = nanRows(z,2);
    s = simAllign.sCoord(inSect,~delCol);

    %replace any missing s-values with mean of that node.
    [mMiss, nMiss] = find(isnan(s));
    for i = 1:length(mMiss)
        s(mMiss(i),nMiss(i)) = nanmean(s(mMiss(i),:));
    end

    %remove mean from each node.
    mz = nanmean(z,2);
    zresid = z - mz;

    grp = groups(~delCol);
    [U,S,V,opts.iSV{r},opts.iSV_orbitVec{r},St] = parallelAnalysis(zresid,100,grp,0.05);

    z2 = SVRecomp(U,S,V,opts.iSV{r});
    
    z2 = z2 + mz;
  
    % join section data for later comparison  
    z2All(inSect,~delCol) = z2;
    sAll(inSect,~delCol) = s;
    
    opts.avgSV = opts.avgSV + size(z,1).*numel(opts.iSV{r});
    
%     find(cumsum(diag(S)) ./ sum(diag(S)) > 0.99,1)
end
missingRows = sum(isnan(z2All),2) == size(z2All,2);
opts.avgSV = opts.avgSV ./ (length(z2All) - sum(missingRows));
 
%constrain profiles
for i = 1:size(z2All,2)
    z2All(:,i) = slopeConstrain(z2All(:,i),opts.maxDiff);
end

skm = nanmean(sAll,2)/1000;


simStats = nodeStats(skm,zAll,truthAllign.(zField));
svdStats = nodeStats(skm,z2All,truthAllign.(zField));
%% 

%Renato's gaussian smoothing.
RL = 10;
sigma = RL/5;
for i = 1:size(zAll,2)
    [smoothProfs(:,i),~] = GaussianAveraging(skm,zAll(:,i), ... 
        simAllign.nWidth(:,i),RL,sigma);
end

smoothStats.(river) = nodeStats(skm,smoothProfs,truthAllign.(zField));

SVDStats.(river) = svdStats;
SVDStats.(river).x = simAllign.easting;
SVDStats.(river).y = simAllign.northing;
SVDStats.(river).meanProf = slopeConstrain(z2All,opts.maxDiff);

SIMStats.(river) = simStats;
SIMStats.(river).x = simAllign.easting;
SIMStats.(river).y = simAllign.northing;

rOpts.(river) = opts;

%% 
%--------------------------------------------------------------------------
% plots
%--------------------------------------------------------------------------
set(0,'defaultAxesFontSize',12,'DefaultAxesFontName','Times New Roman')

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

%original and approx profiles
handle = figure();

ax = tight_subplot(3,1,[.05 .05],[.075 .05],[.1 .05]);

axes(ax(1))
plot(skm,truthAllign.(zField),'k')
ylabel('Elevation (m)')
title('Hydrodynamic Model')
set(ax(1),'XTickLabel',[])

axes(ax(2))
plot(skm,zAll,'k')
ylabel('Elevation (m)')
title('Simulated SWOT')
set(ax(2),'XTickLabel',[])

axes(ax(3))
plot(skm,z2All,'k')
xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Constrained Low-Rank Approximation')

linkaxes(ax);
set(gcf,'Units','normalized','Position',[0.013672 0.013194 0.59922 0.91528])
% set(gcf,'Units','centimeter','Position',[27.269 21.625 18 18])
% set(gca,'XLim',[17 27])
% set(gca,'YLim',[28 35])
% pdfExport(gcf,'/Users/Ted/Documents/DEMPaper/figuresDraft/methods/methods')

figure()
subplot(1,2,1)
scatter(simStats.dailyMAE.*100, svdStats.dailyMAE.*100,'filled')
scatter1to1(gca,'origin');
xlabel('Simulated MAE(cm)')
ylabel('LRA MAE(cm)')
title('MAE')

svdStats.nodePctChange = ((svdStats.totMAE - simStats.totMAE) ./ simStats.totMAE) * 100;
subplot(1,2,2)
text('Units','normalized','position',[0.05 0.7],'String', ... 
    ['Sim: ' num2str(simStats.totMAE.*100) ' cm'], 'FontSize',24)
text('Units','normalized','position',[0.05 0.5],'String', ... 
    ['LRA: ' num2str(svdStats.totMAE.*100) ' cm'], 'FontSize',24)
text('Units','normalized','position',[0.05 0.3],'String', ... 
    ['% change: ' num2str(svdStats.nodePctChange)], 'FontSize',24)
set(gca,'visible','off')
SVDStats.(river).nodePctChange = svdStats.nodePctChange;

SIMStats.(river).totRMSE
SVDStats.(river).totRMSE
((svdStats.totRMSE - simStats.totRMSE) / simStats.totRMSE) * 100
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
% fileName = '/Users/Ted/Documents/GHL_meetings/DAWG_blog/18.06.18/animation.gif';
% delete(fileName)
% figure()
% box on
% set(gcf,'color','white')
% set(gcf,'Units','centimeters','Position',[27.693 28.011 18 6.2442])
% 
% plot(nanmean(s,2)/1000,SVRecomp(U,S,V,1)+mz)
% text('Units','normalized','position',[0.6 0.8],'String','rank: 1', 'FontSize',24)
% xlim = [17 27]; ylim = [28 35];
% set(gca,'XLim',xlim, 'YLim',ylim)
% xlabel('Flow Distance (km)')
% ylabel('Elevation (m)')
% 
% gif(fileName,'DelayTime',0.5,'frame',gcf)
% 
% 
% for i = 2:16
% plot(nanmean(s,2)/1000,SVRecomp(U,S,V,1:i)+mz)
% text('Units','normalized','position',[0.6 0.8],'String',['rank: ' num2str(i)], 'FontSize',24)
% set(gca,'XLim',xlim, 'YLim',ylim)
% xlabel('Flow Distance (km)')
% ylabel('Elevation (m)')
% 
% gif
% end

%% plot vector 2 orbit diff
% oneEV = SVRecomp(U,S,V,1)+mz;
% twoEV = SVRecomp(U,S,V,1:2)+mz;
% 
% figure
% text('Units','normalized','position',[0.6 0.8],'String','rank: 1', 'FontSize',24)
% 
% 
% 
% subplot(2,1,1)
% hold on
% plot(nanmean(s,2)/1000,oneEV(:,grp==2),'k','Linewidth',1)
% plot(nanmean(s,2)/1000,oneEV(:,grp==1),'r','Linewidth',1)
% 
% xlim = [16 26]; ylim = [28 35];
% set(gca,'XLim',xlim, 'YLim',ylim)
% 
% 
% subplot(2,1,2)
% hold on
% plot(nanmean(s,2)/1000,twoEV(:,grp==2),'k','Linewidth',1)
% plot(nanmean(s,2)/1000,twoEV(:,grp==1),'r','Linewidth',1)
% 
% 
% 
% 
% xlim = [16 26]; ylim = [28 35];
% set(gca,'XLim',xlim, 'YLim',ylim)
% 
% 




