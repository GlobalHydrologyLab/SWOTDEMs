clear
close all

% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/PA.mat')
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/PA_mD0.mat')
load('./sacTest.mat')
% load('/Users/Ted/Documents/MATLAB/SWOTDEMs/smoothingTest/PA_mD99.mat')

river = 'Sacramento';
nSlopes = 50;
sigmaMultip = 1;

colors = brewermap(3,'*set1') + 0.1; %brighten
colors(2,2) = colors(2,2) - 0.2; %less green

set(0,'defaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman')

trueAvg = slopeConstrain(SIMStats.(river).z - SIMStats.(river).zErr, ... 
    rOpts.(river).maxDiff);
simAvg = SVDStats.(river).meanProf;
% simAvg = slopeConstrain(smoothStats.(river).z, rOpts.(river).maxDiff);
s = SIMStats.(river).s;

[Rdt,Gd,sId] = relativeSteepness(s.*1000,trueAvg,nSlopes);
trueCutoff = nanmean(Rdt) + sigmaMultip*nanstd(Rdt);
% trueCutoff = nanmedian(Rd) + sigmaMultip*iqr(Rd);
steepMaskTrue = Rdt>trueCutoff;
steepIdTrue = sId(steepMaskTrue);

close all
figure
% ax = tight_subplot(3,1,[.075 .075],[.1 .05],[.1 .1]);

% subplot(2,1,1)
% axes(ax(1))
yyaxis left
set(gca,'YColor',[0 0 0])
plot(s,trueAvg,'k-','Linewidth',1)
hold on
h(1) = plot(s(steepIdTrue),trueAvg(steepIdTrue),'.', 'MarkerSize',10);
ylabel('Elevation (m)')
set(gca,'YLim',[5 40])

yyaxis right
plot(s(sId),Rdt,'Linewidth',1)
hold on
plot(get(gca,'XLim'),[trueCutoff, trueCutoff],'k--')
scaleAxis(gca,1.5)
set(gca,'XTick',0:15:150)
xlabel('Distance (km)')
ylabel('Relative Steepness (m^{-1})')
title('Riffle Pool Classification','FontSize',16)

[Rd,Gd,sId] = relativeSteepness(s.*1000,simAvg,nSlopes);
simCutoff = nanmean(Rd) + sigmaMultip*nanstd(Rd);
% simCutoff = nanmedian(Rd) + sigmaMultip*iqr(Rd);
steepMaskSim = Rd>simCutoff;
steepIdSim = sId(steepMaskSim);

% subplot(3,1,2)
% axes(ax(2))
% yyaxis left
% plot(s,simAvg,'Linewidth',1)
% hold on
% plot(s(steepIdSim),simAvg(steepIdSim),'r.','MarkerSize',10)   
% ylabel('Elevation (m)')
% 
% yyaxis right
% plot(s(sId),Rd,'Linewidth',1)
% hold on
% plot(get(gca,'XLim'),[simCutoff, simCutoff],'k--')
% scaleAxis(gca,1.5)
% set(gca,'XTick',0:15:150)
% % xlabel('Distance (km)')
% ylabel('Relative Steepness (m^{-1})')
% title('SWOT sim')
% set(gcf,'Position',[555 181 1564 1131])


correct = steepMaskTrue + steepMaskSim == 2;
falseNeg = steepMaskTrue - steepMaskSim == 1;
falsePos = steepMaskTrue - steepMaskSim == -1;

% 
% axes(ax(3))
% plot(s,simAvg,'k','Linewidth',1.2)
% hold on
% plot(s(sId(correct)),simAvg(sId(correct)),'.','MarkerSize',10,'Color',colors(1,:))
% plot(s(sId(falsePos)),simAvg(sId(falsePos)),'.','MarkerSize',10,'Color',colors(2,:))
% plot(s(sId(falseNeg)),simAvg(sId(falseNeg)),'.','MarkerSize',10,'Color',colors(3,:))

set(gcf,'Position',[786 745 851 343])
yyaxis left
% set(gca,'YLim',[0 40])
lastClass = 'none';
for i = 1:length(correct)
    sIdx = i + nSlopes;
    
    if correct(i)
        iClass = 'correct';
    elseif falsePos(i)
        iClass = 'falsePos';
    elseif falseNeg(i)
        iClass = 'falseNeg';
    else
        iClass = 'none';
    end
    
    if ~strcmp(iClass,lastClass)
        switch lastClass
            case 'none'
                lastClass = iClass;
                left = (s(sIdx-1) - s(sIdx))/2 + s(sIdx);
                continue                
            case 'correct'              
                plotColor = colors(1,:);               
            case 'falsePos'                
                plotColor = colors(3,:);     
            case 'falseNeg'
                plotColor = colors(2,:);
        end
        
        right = (s(sIdx+1) - s(sIdx))/2 + s(sIdx);
        rectangle('Position',[left, 5, right-left, 5], ... 
                    'FaceColor',plotColor, 'EdgeColor','none')
        
        left = (s(sIdx-1) - s(sIdx))/2 + s(sIdx); 
        lastClass = iClass;
    end
    
    
    
%     if correct(i)
%         rectangle('Position',[left, 5, right-left, 5], ... 
%                'FaceColor',colors(1,:), 'EdgeColor','none')
%     end
%     if falsePos(i)
%         rectangle('Position',[left, 5, right-left, 5], ... 
%                'FaceColor',colors(3,:), 'EdgeColor','none')
%     end
%     if falseNeg(i)
%         rectangle('Position',[left, 5, right-left, 5], ... 
%                'FaceColor',colors(2,:), 'EdgeColor','none')
%     end
end

%klugey fix to get rectangles in legend.
h(end+1) = plot(0,0,'s', 'MarkerSize',10,'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
h(end+1) = plot(0,0,'s', 'MarkerSize',10,'Color',colors(2,:),'MarkerFaceColor',colors(2,:));
h(end+1) = plot(0,0,'s', 'MarkerSize',10,'Color',colors(3,:),'MarkerFaceColor',colors(3,:));

legend(h,{'True steep points','Correct','False negative','False positive'})
set(gca, 'Layer','top')

% set(gcf,'Position',[1054 530 672 477])
% legend('Profile','Correct','False Positive','False Negative')
% set(gca,'XTick',0:15:150)
% set(gca,'YLim',[5 40])
% 
% title('Evaluation','FontSize',16)
% ylabel('Elevation (m)')
% xlabel('Distance (km)')

fprintf(['correct:\t' num2str(sum(correct)) '\n'])
fprintf(['false negative:\t' num2str(sum(falseNeg)) '\n'])
fprintf(['false positive:\t' num2str(sum(falsePos)) '\n'])

% pdfExport(gcf,'/Users/Ted/Documents/DEMPaper/figuresDraft/rifflePool/rifflePoolBar')



