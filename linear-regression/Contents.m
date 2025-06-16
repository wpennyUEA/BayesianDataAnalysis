%
% BAYESIAN LINEAR REGRESSION
%
% Bayesian Inference with Multivariate Gaussian Distributions
% Bayesian Inference with Multivariate T Distributions
% Savage-Dickey Bayes Factors 
% Implements Default Bayes Factors (mostly Schwarzkopf Code)
%
% ---------------------------------------------------------------------
% bayes_glm_ttest1.m                One-sample t-tests
% bayes_glm_anova1.m                One-way ANOVA
% bayes_glm_regression.m            Regression
%
% ---------------------------------------------------------------------
% demo_regression.m                 SD-MVT versus Default
% demo_anova1_chick_weights.m       Compare to R code
% demo_anova1_compare.m             SD-MVT versus SD-Gauss versus Default [3,4]
% demo_anova1_confounds.m           ANOVA's with confounds
% demo_ttest_standard_scale.m       Compare SDT versus Default
% demo_ttest_arbitrary_scale.m      Compare SDT versus Default
%
% ---------------------------------------------------------------------
% bayes_linear_estimate.m           Generic Bayesian inference for Linear Model
% bayes_linear_test.m               Savage-Dickey (SD) Bayes Factors
%
% bayes_tlinear_estimate.m          Multivariate-T (MVT) Distributions [1,2]
% bayes_tlinear_test.m              SD with Multivariate-T (MVT) Distributions
%
% bayes_glinear_estimate.m          Gaussian Distributions
% bayes_glinear_test.m              SD with Gaussian Distributions
%
% ---------------------------------------------------------------------
% linregbf.m                        Default Bayes Factors for Linear Regression [5]
% linregbf_test.m                   Validates linregbf.m against corrbf.m for single regressor 
%
% ---------------------------------------------------------------------
% glm_anova1.m                      Classical Inference with SPM code
% glm_test_hypothesis.m             SPM contrast matrix approach
% make_anova_design.m               Create design matrices
%
% ---------------------------------------------------------------------
% REFERENCES
%
% [1] Bernardo and Smith, Bayesian Theory, Wiley, 2000. See Appendix A.
%
% [2] W. Penny. Bayesian General Linear Models with T-Priors. 
% Technical Report, UEA, 2020. 
% 
% [3] Rouder et al., 2009 Bayesian t tests for accepting and rejecting the null
% hypothesis. Psych Bull Review, 16(2),225-237.
%
% [4] Wetzels et al, 2012. A default Bayesian hypothesis test for ANOVA
% designs. The American Statistician, 66(2),104-111.
%
% [5] Rouder and Morey, Default Bayes Factors for Model Selection 
% in Regression. Multivariate Behavioural Research, 2013.
%
