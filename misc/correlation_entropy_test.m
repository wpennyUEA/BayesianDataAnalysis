
clear all
close all

d=10;
r=[0:0.1:0.9];
N=length(r);
for n=1:N,
    R=r(n)*ones(d,d);
    for i=1:d,
        R(i,i)=1;
    end
    H(n)=gaussian_entropy(R);
end

figure
plot(r,H);
grid on
xlabel('Pairwise Correlation');
ylabel('Entropy of Correlation Matrix');
