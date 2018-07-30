clear
close all

% load('Sacramento/DEMProfiles/dems.mat')
% load('./smoothingTest/PA_geoidFix.mat')

load('Figures/DEMCompare/interpConstrainDEMs.mat')

nSlopes = 50;
sigmaMultip = 1;
maxDiff = 0;

colors = brewermap(3,'*set1') + 0.1; %brighten
colors(2,2) = colors(2,2) - 0.2; %less green

set(0,'defaultAxesFontSize',12,'DefaultAxesFontName','Times New Roman')

trueAvg = slopeConstrain(Sac.avgTruth, maxDiff);
simAvg = Sac.SWOT;
s = Sac.skm;

[Rdt,Gdt,sId] = relativeSteepness(s.*1000,trueAvg,nSlopes);
trueCutoff = nanmean(Rdt) + nanstd(Rdt);

steepMaskTrue = Rdt>trueCutoff;
steepIdTrue = sId(steepMaskTrue);

DEMs = {'SWOT','SRTM','MERIT','ASTER','NED','lidar'};

for i = 1:numel(DEMs)
    comp.(DEMs{i}) = relativeSteepnessComp(Sac.skm,Sac.(DEMs{i}),nSlopes,steepMaskTrue);
end

DEMs = fields(comp);

for i = 1 : numel(DEMs)
    stats(i,:) = comp.(DEMs{i}).Stats;
end

Correct = stats(:,1);
falseNegative = stats(:,2);
falsePositive = stats(:,3);
errTable = table(Correct,falseNegative,falsePositive,'RowNames',DEMs)
%%

close all
figure

yyaxis left
set(gca,'YColor',[0 0 0])
plot(s,trueAvg,'k-','Linewidth',1)
hold on
h(1) = plot(s(steepIdTrue),trueAvg(steepIdTrue),'.', 'MarkerSize',10);
ylabel('Elevation (m)')
set(gca,'YLim',[0 40])

yyaxis right
plot(s(sId),Rdt,'Linewidth',1)
hold on
plot(get(gca,'XLim'),[trueCutoff, trueCutoff],'k--')
scaleAxis(gca,1.5)
set(gca,'XTick',0:15:150)
xlabel('Distance (km)')
ylabel('Relative Steepness (m^{-1})')
title('Steepness Classification')


correct = comp.SWOT.correct;
falseNeg = comp.SWOT.falseNeg;
falsePos = comp.SWOT.falsePos;

set(gcf,'Units','centimeters','Position',[27.693 26.247 18 8])
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
        rectangle('Position',[left, 0, right-left, 5], ... 
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
h(end+1) = plot(-5,-5,'s', 'MarkerSize',10,'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
h(end+1) = plot(-5,-5,'s', 'MarkerSize',10,'Color',colors(2,:),'MarkerFaceColor',colors(2,:));
h(end+1) = plot(-5,-5,'s', 'MarkerSize',10,'Color',colors(3,:),'MarkerFaceColor',colors(3,:));

legend(h,{'True steep points','Correct','False negative','False positive'})
set(gca, 'Layer','top')

% set(gcf,'Position',[1054 530 672 477])
% legend('Profile','Correct','False Positive','False Negative')
% set(gca,'XTick',0:15:150)
% set(gca,'YLim',[5 40])

% 
% fprintf(['correct:\t' num2str(sum(correct)) '\n'])
% fprintf(['false negative:\t' num2str(sum(falseNeg)) '\n'])
% fprintf(['false positive:\t' num2str(sum(falsePos)) '\n'])

% pdfExport(gcf,'/Users/Ted/Documents/DEMPaper/figuresDraft/rifflePool/rifflePoolBar')



