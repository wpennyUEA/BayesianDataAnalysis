function [F] = blr_compare (model,name,t,opt)
% Compare sets of regressors using model evidence
% FORMAT [F] = blr_compare (model,name,t,opt)
%
% model     (m).x       design matrix for model m
% name      {m}         name of model m
% t         labels
% opt       optimisation params
%
% F         (m) log evidence of model m

M = length(model);
for m=1:M,
    if isempty(model(m).x)
        X = [];
    else
        X = zmuv(model(m).x); % Ensure all inputs have zero mean and unit variance
    end
    M = blr_fit(X,t,opt);
    F(m) = M.F;
end

figure
bar(F-min(F));
grid on
set(gca,'XTickLabel',name);
xlabel('Model');
ylabel('Log Evidence');
