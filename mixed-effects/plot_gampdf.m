

mean = 5;
var = 0.3;

% Jean's parameterisation
b = mean/var;
a = b*mean;

% Parameterisation for spm_kl_gamma
my_b = 1/b;
c = a;

x=[0.1:0.1:10];

[p,m,v] = mygampdf(x,my_b,c);

figure
plot(x,p);
title(sprintf('Mean = %1.3f Var = %1.3f',m,v));