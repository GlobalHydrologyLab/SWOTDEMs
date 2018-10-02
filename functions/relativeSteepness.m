function [Rd,Gd,xId] = relativeSteepness(x,z,nSlopes)
%calculating relative steepness from Hayakawa and Oguchi 2006

nx = numel(z);

minidx = 1+nSlopes;
maxidx = nx-nSlopes;

[leftA,leftB] = meshgrid(flip(1:nSlopes),0:(maxidx-minidx));
left = leftA + leftB;

right = fliplr(left)+minidx;

d = x(left)-x(right);
Gd = -1.*(z(left)-z(right)) ./ d;

for i = 1:size(Gd,1)
    Rd(i,:) = polyfit(d(i,:),Gd(i,:),1);
end

Rd(:,2) = [];
xId = minidx:maxidx;

end