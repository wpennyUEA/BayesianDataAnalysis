function [R] = svm_classify(X,Y,V,perm)
% V-fold Cross Validation using a Support Vector Machine
% FORMAT [R] = svm_classify(X,Y,V,perm)
%
% INPUTS:
% X     [N x d] data matrix
% Y     [N x 1] element string vector of category labels (K different categories)
% V     Number of folds, default = 10
% perm  1 to randomly permute samples before cross-validation (default 0)
%
% OUTPUTS:
% R     .CM         confusion matrix
%       .order      of variables in CM
%       .corr_rate  correct rate
%       .err_rate   error rate
%
%       For K=2:
%       .m_auc      Area Under Curve (AUC)
%       .sem_auc    AUC Standard Error 
%
% Requires Matlab's statistics toolbox


try V=V; catch V = 10; end
try perm=perm; catch perm=0; end

% Add "bias input"
% [N,d] = size(X);
% X = [X,ones(N,1)];

Y = Y(:);
N = size(X,1);
labels = unique(Y);
K = length(labels);

if K==2, compute_auc = 1; else compute_auc=0; end
if perm, ind = randperm(N); else ind = 1:N; end

block_size = floor(N/V);
cm_tot = zeros(K,K);
for i=1:V,
    start = (i-1)*block_size+1;
    stop = start + block_size-1;
    tr = ind;
    tr(start:stop)=[];
    te = ind(start:stop);
    
    if K==2
        %S = fitcsvm(X(tr,:),Y(tr));
        S = fitclinear(X(tr,:),Y(tr));
        %S = fitclinear(X(tr,:),Y(tr),'Regularization','lasso');
        %S = fitclinear(X(tr,:),Y(tr),'Regularization','ridge');
    else
        S = fitcecoc(X(tr,:),Y(tr));
        %S = fitcecoc(X(tr,:),Y(tr),'Coding','onevsall');
    end

    [label,score] = predict(S,X(te,:));
    if K==2
        %tmp=strcmp(labels{1},Y(te));
        tmp=strcmp(label(:),Y(te));
        i0=find(tmp==1);
        i1=find(tmp==0);
        if isempty(i0)
            auc_svm(i)=0.5;
        elseif isempty(i1)
            auc_svm(i)=1;
        else
            auc_svm(i) = auc_ranksum(score(i0,2),score(i1,2));
        end
    end
    %[cm,order] = confusionmat(Y(te),label);

    [cm,order] = confusionmat(Y(te),label,'order',S.ClassNames);

    cm_tot = cm_tot + cm;

end

R = [];
if K==2
    R.m_auc = mean(auc_svm);
    R.sem_auc = std(auc_svm)/sqrt(V);
end
R.CM = cm_tot;
R.order = S.ClassNames;
R.corr_rate = sum(diag(R.CM))/sum(sum(R.CM));
R.err_rate = 1 - R.corr_rate;


