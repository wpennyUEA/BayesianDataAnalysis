


% Jean's parameterisation
b = 1e-4;
a = 1e-3;

m = a/b
v = a/(b^2)

% Parameterisation for spm_kl_gamma
my_b = 1/b;
c = a;

x=[0.1:0.1:100];

[p,m,v] = mygampdf(x,my_b,c);

figure
plot(x,p);
title(sprintf('Mean = %1.3f Var = %1.3f',m,v));