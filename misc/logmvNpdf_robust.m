function [logp] = logmvNpdf_robust(z,Mu,V)
% Log pdf of multivariate Normal distribution
% FORMAT [logp] = logmvNpdf_robust(z,Mu,V)
%
% z  - ordinates
% Mu - mean (a d-vector)
% V  - d x d variance-covariance matrix
%

%-Condition arguments
%--------------------------------------------------------------------------
if nargin<1,   logp=[]; return, end
if isempty(z), logp=[]; return, end
if nargin<2,   Mu=0;           end

%-Check Mu, make a column vector, get dimension
%--------------------------------------------------------------------------
if min(size(Mu)) > 1, error('Mu must be a vector'); end
Mu = Mu(:)';
d  = length(Mu);

if nargin<3, V=eye(d); end

%-Size & range checks
%--------------------------------------------------------------------------
if any(any(V~=V')),     error('V must be symmetric'); end
if any(size(V)~=[d,d]), error('V wrong dimension');   end

%-Computation
%--------------------------------------------------------------------------
if d==1
    %-Simpler computation for univariate normal
    %----------------------------------------------------------------------
    logp = -(z - Mu).^2/(2*V)-log(sqrt(2*pi*V));
else
    if size(z,1) ~= d, error('z wrong dimension'), end
    z   = z - Mu(:)*ones(1,size(z,2));
    logp = -0.5*d*log(2*pi) - 0.5*spm_logdet(V) - 0.5*z'*inv(V)*z;
    %logp = -0.5*sum((sqrtm(inv(V))*z).^2)-log((2*pi)^(d/2)*sqrt(det(V)));
end

