clear
close all

load('SWOTDEMs/Figures/DEMCompare/poInterp.mat')

smoothWindow = 51;
cWindow = 101; % 20km window

truth = calcConvexity(Po.skm,Po.avgTruth,smoothWindow,cWindow);
swot = calcConvexity(Po.skm,Po.SWOT,smoothWindow,cWindow);
srtm = calcConvexity(Po.skm,Po.SRTM,smoothWindow,cWindow);
merit = calcConvexity(Po.skm,Po.MERIT,smoothWindow,cWindow);
aster = calcConvexity(Po.skm,Po.ASTER,smoothWindow,cWindow);
tinitaly = calcConvexity(Po.skm,Po.TINITALY,smoothWindow,cWindow);


sc = 0.3;

figure
hold on
box on

xl = [0,130];
set(gca,'XLim',xl)


set(gca,'YLim',[-sc/2,5.5*sc])
set(gca,'YTick',-sc/4:sc/2:5.25*sc)

plot(xl,[(0:5)*sc; (0:5)*sc],'k--')
plot(xl,[((sc/2):sc:(5*sc)); ((sc/2):sc:(5*sc))],'k-')

plot(Po.skm,truth,'k','Linewidth',2)
plot(Po.skm,swot + 1*sc,'Linewidth',2)
plot(Po.skm,srtm + 2*sc,'Linewidth',2)
plot(Po.skm,merit + 3*sc,'Linewidth',2)
plot(Po.skm,aster + 4*sc,'Linewidth',2)
plot(Po.skm,tinitaly + 5*sc,'Linewidth',2)

xlabel('Flow Distance (m)')

set(gca,'YTickLabel',reshape((ones(5,1)*[-sc/2,sc/2])',[10,1]))

yyaxis right
set(gca,'YTick',1/12 : 1/6 : 11/12)
set(gca,'YTickLabel',{'Hydrodynamic Model','SWOT','SRTM','MERIT','ASTER','TINITALY'})
set(gca,'YColor',[0 0 0])


set(gcf,'Units','centimeters','Position',[28 19.121 17.992 15.134])
pdfExport(gcf,'/Users/Ted/Documents/DEMPaper/figuresDraft/poConvexities/demCurvatures_10km')

function [c] = calcConvexity(s,z,smoothWindow,cWindow)
    n = length(s);
    z = smooth(z,smoothWindow);
    
    c = nan(n,1);
    for i = 1+((cWindow-1)/2) : n-((cWindow-1)/2)-1
        span = [i-((cWindow-1)/2), i, i+((cWindow-1)/2)];
        c(i) = diff(z(span),2)/range(s(span));
    end

end