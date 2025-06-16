
clear all
close all

p=6;
timepoints=512;
channels=64;
X=randn(timepoints,channels);
tic;
mar=spm_mar(X,p);
toc

%--- 98 seconds on Will's desktop 8G RAM
% with
% p=6;
% timepoints=512;
% channels=32;

%--- 28 seconds on Will's newer desktop 8G RAM
% with
% p=6;
% timepoints=512;
% channels=32;

%--Didn't finish even overnight !!! on Will's newer desktop 8G RAM
% with
% p=6;
% timepoints=512;
% channels=64;