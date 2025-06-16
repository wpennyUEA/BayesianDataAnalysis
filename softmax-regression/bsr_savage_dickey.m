function [logBF,z] = bsr_savage_dickey (bsr)
% Computer Bayes factors in favour of pruning parameters
% FORMAT [logBF,z] = bsr_savage_dickey (bsr)
%
% bsr      data structure
%
% logBF    logBF(k,d) Bayes Factor in favour of hypothesis 
%          that w(k,d)is non-zero
% z        z(k,d) is z-score from posterior distribution

K = length(bsr.class);
D = length(bsr.class(1).w);
for k=1:K,
    for d=1:D-1,
        wkd0 = bsr.class(k).w0(d);
        wkd = bsr.class(k).w(d);
        rkd0 = bsr.class(k).R0(d,d);
        rkd = bsr.class(k).R(d,d);
        prior = spm_Npdf(0,wkd0,1/rkd0);
        post = spm_Npdf(0,wkd,1/rkd);
        logBF(k,d) = log(prior/post);
        z(k,d) = abs(wkd)/sqrt(1/rkd);
    end
end