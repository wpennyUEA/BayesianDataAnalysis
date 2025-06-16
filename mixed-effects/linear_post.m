function [Ep,Cp,L] = linear_post (M,U,Y)
% Analytic posterior for linear regression
% FORMAT [Ep,Cp,L] = linear_post (M,U,Y)
% 
% M     Model Structure
% U     Inputs
% Y     Data
%
% Ep    Posterior mean
% Cp    Posterior covariance
% L     Log evidence

if isstruct(Y)
    Y=Y.y;
end

% Posteriors
ipC=inv(M.pC);
iCe=inv(M.Ce);
X=U.X;
iCp=X'*iCe*X+ipC;
Cp=inv(iCp);
Ep=Cp*(X'*iCe*Y+ipC*M.pE);

% Log evidence
yhat=X*Ep;
T=length(yhat);
ey=Y-yhat;

L = -0.5*T*spm_logdet(M.Ce) - 0.5*T*log(2*pi);
L = L - 0.5*trace(ey'*iCe*ey);

ew = M.pE-Ep;
L = L - 0.5*ew'*ipC*ew - 0.5*spm_logdet(M.pC) + 0.5*spm_logdet(Cp);
