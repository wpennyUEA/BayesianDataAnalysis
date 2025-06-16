clear all
close all

disp('Two-dimensional data with three clusters');
disp('Incremental GMM learning (fully online)');

% Use data from file
load yrep
figure
plot(y(:,1),y(:,2),'x');
hold on
   
x = y';
s0 = mean(std(y));

mix = igmm_init (x(:,1),s0);
[d,N]=size(x);
for n=2:N,
    %mix = igmm_update_diag (mix,x(:,n));
    mix = igmm_update (mix,x(:,n));
end

mix.m = mix.M;

for i=1:mix.M,
   mix.state(i).C = inv(mix.state(i).Lambda);
   plot(mix.state(i).m(1),mix.state(i).m(2),'rx');
end
hold on
spm_mix_plot2d(mix,[-2 10 -2 10],1,'r',0.4,0.5);
set(gca,'FontSize',18);
grid on


