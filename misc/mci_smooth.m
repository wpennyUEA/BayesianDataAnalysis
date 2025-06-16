function [y,std_dev] = mci_smooth (x, win)
% FORMAT [y,std_dev] = mci_smooth (x, win)
% Moving Average filter without edge effects
%
% x     time series to be smoothed
% win   window size to average over (number of samples)
%
% y         smoothed time series
% std_dev   std deviation thereof (due to smoothing)
%
% y(t) is the mean value over left(t) and right(t) samples
% where left(t)=0 up to win/2 and right(t)=0 up to win/2
% and a variable size of left and right windows is use to avoid
% edge effects.

T=length(x);
s=[1:T];
for t=1:T,
    left=find(s>(t-win/2) & s<t);
    right=find(s>t & s < (t+win/2));
    
    ind=[left,t,right];
    Nind=length(ind);
    y(t)=mean(x(ind));
    std_dev(t) = std(x(ind))/sqrt(Nind);
end