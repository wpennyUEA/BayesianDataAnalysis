function [y] = norm_space (x,msd)
% Project data from original space to the zmuv space
% FORMAT [y] = norm_space (x,msd)

I = size(x,1);
for i=1:I,
    y(i,:) = (x(i,:)-msd(i,1))/msd(i,2);
end