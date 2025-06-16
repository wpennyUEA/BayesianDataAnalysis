function [x] = org_space (y,msd)
% Project data from zmuv space to original space
% FORMAT [x] = org_space (y,msd)

I = size(y,1);
for i=1:I,
    x(i,:) = msd(i,2)*y(i,:)+msd(i,1);
end