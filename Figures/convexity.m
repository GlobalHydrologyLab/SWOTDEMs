clear
close all
addpath('/Users/Ted/Documents/MATLAB/cited_functions/fitcircle')

load('./smoothingTest/PA_geoidFix.mat')
set(0,'defaultAxesFontSize',12,'DefaultAxesFontName','Times New Roman')


md = 0.01;


s = SVDStats.Po.s;
z = slopeConstrain(SVDStats.Po.z,md);
zt = slopeConstrain(SIMStats.Po.z - SIMStats.Po.zErr,md);
ztp = SIMStats.Po.z - SIMStats.Po.zErr;


n = length(s);
smoothWindow = 51;
cWindow = 101; % 20km window

z = smooth(z,smoothWindow);
% z = GaussianAveraging(s,z,z,20,20/5);

Z = NaN(n,2);
R = NaN(n,1);
c = NaN(n,1);

for i = 1+((cWindow-1)/2) : n-((cWindow-1)/2)-1
%     span = i-((cWindow-1)/2) : i +((cWindow-1)/2)+1;
%     [Z(i,:),R(i)] = fitcircle([s(span) z(span)]);
    span = [i-((cWindow-1)/2), i, i+((cWindow-1)/2)];
    c(i) = diff(z(span),2)/range(s(span));
end

%% 


% figure
% yyaxis left
% plot(s,z)
% ylabel('Elevation (m)')
% yyaxis right
% plot(s,c)
% ylabel('Convexity')

colors = brewermap(100,'RdBu');
cNorm = c ./ max(abs(c)) /2;

cScale = round(cNorm*100)+50;

figure
hold on
box on

for i = 1+((cWindow-1)/2) : n-((cWindow-1)/2)-1
    x1 = (s(i-1)+s(i)) /2;
    y1 = (z(i-1)+z(i)) /2;
    x2 = (s(i)+s(i+1)) /2;
    y2 = (z(i)+z(i+1)) /2;
    plot([x1 x2],[y1 y2], 'Linewidth', 4, ... 
        'Color',colors(cScale(i),:))
end
colormap(colors);
h = colorbar;

cx = [0 0.2 0.4 0.6 0.8 1];
cMax = max(abs(c));
cy = round([-cMax -cMax/(3/5) -cMax/(1/5) cMax/(1/5) cMax/(3/5) cMax],3);

h.Ticks = cx;
h.TickLabels = cy;
ylabel(h,'Concavity (m/km^2)')

xlabel('Flow Distance (km)')
ylabel('Elevation (m)')
title('Concavity of the Po River Profile')
set(gca,'YLim',[2.5 17.5],'XLim',[2.5,130])


set(gcf,'Units','centimeters','Position',[27.693 27.093 18 7.1614])
% pdfExport(gcf,'/Users/Ted/Documents/DEMPaper/figuresDraft/poConvexities/poCurvature_10km')









