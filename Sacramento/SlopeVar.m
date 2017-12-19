%looking at ways to show relationship between slope and width.

clear
close all
load('Sacramento/SlopeVarData.mat')

k=21; %window size

st = truthAvg.sCoord;
ss = simAvg.sCoord;

zt = truthAvg.geoHeight;
zs = simAvg.geoHeight;

xt = truthAvg.easting;
xs = simAvg.easting;

yt = truthAvg.northing;
ys = simAvg.northing;

w = simAvg.nWidth;

%%
%node averaged plot.
figure(1)
plot(movstd(zt,k),movstd(zs,k),'.')
line([0,0.8],[0,0.8])
xlabel('windowed std. dev. of true elevations (m)')
ylabel('windowed std. dev. of simulated elevations (m)')
% grid on


%plot s,z profiles with color bar for movstd(z) values




%% test detection of meanders
close all

transParam = [1 3 11 1600 200]';
[sn,~,cl_out] = xy2sn([xt,yt],[xt,yt],transParam);

plot(cl_out(:,1),cl_out(:,2))
hold on
plot(xt,yt)
