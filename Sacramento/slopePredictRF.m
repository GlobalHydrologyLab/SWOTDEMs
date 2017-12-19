%--------------------------------------------------------------------------
clear;clc;close all

%--------------------------------------------------------------------------
% Load an example dataset provided with matlab

% load chemical_dataset
% In = chemicalInputs';
% Out = chemicalTargets';

%load my own data. see what can predict slopes best.
k = 7; %39 node window approximates 1km^2 area for average width of Sac.
load('/Users/Ted/Documents/MATLAB/SWOT_DEMs/Sacramento/SlopeVarData.mat')

% % variance of elevations (biased by overall slope)
% Out = movstd([truth.geoHeight],k,0,1,'omitnan');
% Out = reshape(Out,[],1); %vector

% slope of elevations, with an extra appended to keep same dimensions.
Out = diff([truth.geoHeight],1,1)./diff([truth.sCoord],1,1).*-1;
Out(end+1,:) = Out(end,:);
Out = movstd(Out,k,0,1,'omitnan');
Out = reshape(Out,[],1);


% testVar{1} = 'local variance of simulated heights (m)';
% tmp = movstd([simulated.geoHeight],k,0,1,'omitnan');
% In(:,1) = reshape(tmp,[],1);

testVar{1} = 'local std. dev. of simulated slopes (m)';
tmp = diff([simulated.geoHeight],1,1)./diff([simulated.sCoord],1,1).*-1;
tmp(end+1,:) = tmp(end,:);
tmp = movstd(tmp,k,0,1,'omitnan');
In(:,1) = reshape(tmp,[],1);

testVar{2} = 'local std. dev of simulated widths (m)';
tmp = movstd([simulated.nWidth],k,0,1,'omitnan');
In(:,2) = reshape(tmp,[],1);

testVar{3} = 'local mean of planform angular deflection (deg/m)';
for i = 1:length(simulated)
    tmp = angDeflection(simulated(i).easting,simulated(i).northing);
    aD(:,i) = [tmp(1); tmp; tmp(end)]';
end
aD = smooth(aD, k, 'moving');
In(:,3) = reshape(aD,[],1);

% testVar{4} = 'windowed mean of CL curvature from riverKrige';
% transParam = [1 3 11 717 200]';
% for i = 1:length(simulated)
%     tmp = curvature([simulated(i).easting,simulated(i).northing],transParam);
%     cur(:,i) = tmp(:,2);
% end
% In(:,4) = reshape(cur,[],1);

%--------------------------------------------------------------------------
% Find capabilities of computer so we can best utilize them.

% Find if gpu is present
ngpus=gpuDeviceCount;
disp([num2str(ngpus) ' GPUs found'])
if ngpus>0
    lgpu=1;
    disp('GPU found')
    useGPU='yes';
else
    lgpu=0;
    disp('No GPU found')
    useGPU='no';
end

% Find number of cores
ncores=feature('numCores');
disp([num2str(ncores) ' cores found'])

% Find number of cpus
import java.lang.*;
r=Runtime.getRuntime;
ncpus=r.availableProcessors;
disp([num2str(ncpus) ' cpus found'])

if ncpus>1
    useParallel='yes';
else
    useParallel='no';
end

[archstr,maxsize,endian]=computer;
disp([...
    'This is a ' archstr ...
    ' computer that can have up to ' num2str(maxsize) ...
    ' elements in a matlab array and uses ' endian ...
    ' byte ordering.'...
    ])

% Set up the size of the parallel pool if necessary
npool=ncores;

% Opening parallel pool
if ncpus>1
    tic
    disp('Opening parallel pool')
    
    % first check if there is a current pool
    poolobj=gcp('nocreate');
    
    % If there is no pool create one
    if isempty(poolobj)
        command=['parpool(' num2str(npool) ');'];
        disp(command);
        eval(command);
    else
        poolsize=poolobj.NumWorkers;
        disp(['A pool of ' poolsize ' workers already exists.'])
    end
    
    % Set parallel options
    paroptions = statset('UseParallel',true);
    toc
    
end

%--------------------------------------------------------------------------
tic
leaf=5;
ntrees=200;
fboot=1;
surrogate='on';
disp('Training the tree bagger')
b = TreeBagger(...
        ntrees,...
        In,Out,... 
        'Method','regression',...
        'oobvarimp','on',...
        'surrogate',surrogate,...
        'minleaf',leaf,...
        'FBoot',fboot,...
        'Options',paroptions,...
        'InBagFraction',0.7...
    );
toc

%--------------------------------------------------------------------------
% Estimate Output using tree bagger
disp('Estimate Output using tree bagger')
x=Out;
y=predict(b, In);
name='Bagged Decision Trees Model';
toc

%--------------------------------------------------------------------------
% calculate the training data correlation coefficient
cct=corrcoef(x,y,'Rows','complete');
cct=cct(2,1);

%--------------------------------------------------------------------------
% Create a scatter Diagram
disp('Create a scatter Diagram')

% plot the 1:1 line
plot(x,x,'LineWidth',3);

hold on
scatter(x,y,'filled');
hold off
grid on

set(gca,'FontSize',18)
xlabel('Actual','FontSize',25)
ylabel('Estimated','FontSize',25)
title(['Training Dataset, R^2=' num2str(cct^2,2)],'FontSize',30)

drawnow

fn='ScatterDiagram';
fnpng=[fn,'.png'];
print('-dpng',fnpng);

%--------------------------------------------------------------------------
% Calculate the relative importance of the input variables
tic
disp('Sorting importance into descending order')
weights=b.OOBPermutedVarDeltaError;
[B,iranked] = sort(weights,'descend');
toc

%--------------------------------------------------------------------------
disp('Plotting a horizontal bar graph of sorted labeled weights.') 

%--------------------------------------------------------------------------
figure
barh(weights(iranked));
xlabel('Variable Importance');
ylabel('Variable Rank');
title('Relative importance of SWOT simulated variables in estimating true slope');
hold on

%--------------------------------------------------------------------------
grid on 
xt = get(gca,'XTick');    
xt_spacing=unique(diff(xt));
xt_spacing=xt_spacing(1);    
yt = get(gca,'YTick');    
ylim([0.25 length(weights)+0.75]);
xl=xlim;
xlim([0 2.5*max(weights)]);

%--------------------------------------------------------------------------
% Add text labels to each bar
for ii=1:length(weights)
    text(...
        max([0 weights(iranked(ii))+0.02*max(weights)]),ii,...
        [testVar{iranked(ii)}]);
end

%--------------------------------------------------------------------------
set(gca,'FontSize',16)
set(gca,'XTick',0:2*xt_spacing:1.1*max(xl));
set(gca,'YTick',yt);
set(gca,'TickDir','out');
set(gca, 'ydir', 'reverse' )
set(gca,'LineWidth',2);   
drawnow

%--------------------------------------------------------------------------
fn='RelativeImportanceInputs';
fnpng=[fn,'.png'];
print('-dpng',fnpng);

%--------------------------------------------------------------------------
% Ploting how weights change with variable rank
disp('Ploting out of bag error versus the number of grown trees')

figure
plot(b.oobError,'LineWidth',2);
xlabel('Number of Trees','FontSize',30)
ylabel('Out of Bag Error','FontSize',30)
title('Out of Bag Error','FontSize',30)
set(gca,'FontSize',16)
set(gca,'LineWidth',2);   
grid on
drawnow
fn='EroorAsFunctionOfForestSize';
fnpng=[fn,'.png'];
print('-dpng',fnpng);

