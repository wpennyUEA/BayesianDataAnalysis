
clear all
close all

% Use data from file
load yrep
figure
plot(y(:,1),y(:,2),'x');
hold on
   
m_model=5;
disp('Two-dimensional data with three clusters');
disp(sprintf('Assumed model has %d clusters',m_model));
disp('VB GMM code');

vbmix=spm_mix(y,m_model);

for i=1:m_model,
   plot(vbmix.state(i).m(1),vbmix.state(i).m(2),'rx');
end
hold on
spm_mix_plot2d(vbmix,[-2 12 -2 12],1,'r',0.4,0.5);
set(gca,'FontSize',18);


