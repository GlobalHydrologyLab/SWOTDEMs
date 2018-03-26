function [peaks, iPeaks] = slidePeaks(y, multip, plot)
%slidePeaks.m uses the built in function 'findpeaks.m' but adds ability to
%use peak prominence as a multiple of local mean, instead of static
%value. takes in signal data and variables for window size and mean
%multiplier.

% y = signal
% multip = 0-1, y prominence thresholding as percent of y value
% plot = plot boolean

if ~exist('multip','var')
    multip = 1/3;
end

if ~exist('p','var')
    plot = 0;
end

[Ypk,Xpk,~,Ppk] = findpeaks(y);

peakThresh = (Ypk - Ppk) .* multip;

meetsCrit = Ppk >= peakThresh;


iPeaks = Xpk(meetsCrit);
peaks = Ypk(meetsCrit);


if plot == 1
    plot(y,'Linewidth',2)
    hold on
    plot(iPeaks,peaks,'kv','MarkerFaceColor','k')
    hold off
end
    
end

