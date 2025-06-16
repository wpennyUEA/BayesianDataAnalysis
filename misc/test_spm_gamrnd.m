
clear all
close all

a=2;
b=3;

m=a*b
var=b^2*a

N=1000;
r=spm_gamrnd(a,b,N,1);

mean(r)
std(r)^2