function [vbmix logev mix] = spm_mix_demo1d (data, maxcomps, verbosity)
% Demonstrate use of spm_mix on 1D data
% FORMAT [vbmix logev mixdata] = spm_mix_demo1d (data, maxcomps, plotfits)
%
% data      - either scalar number of clusters to simulate or your own data
% maxcomps  - maximum number of components in mixture model to consider
% verbosity - 0 = silent, 1 = basic output (with figures), 2 = full output
%
% vbmix     - cell array of fitted mixtures for all numbers of components
% logev     - log evidence for each number of components
% mix       - mix structure for simulated mixtures if scalar data given
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Will Penny & Ged Ridgway
% $Id: spm_mix_demo1d.m 3997 2010-07-15 12:38:24Z ged $


clear all
close all

% This script fits GMMs to 1D data using Variational Bayesian inference

% Defaults
if nargin < 3
    verbosity = 1; % Verbosity level
end
if nargin < 2
    maxcomps = 5; % Maximum number of components to try
end
if nargin < 1
    data = 3; % Simulate data with 3 clusters if no input
end

% If scalar, simulate 1D data from Gaussian mixture with data clusters
if isscalar(data)
    mix.m = data;   % Number of clusters
    if verbosity > 0
        fprintf('Simulating data with %d clusters\n', mix.m)
    end
    means = 0:5:5*(mix.m - 1);  % Define means 

    for m = 1:mix.m
        mix.state(m).prior = 1 / mix.m; % Equal prior weights for each component
        mix.state(m).m = means(m); % Set mean
        mix.state(m).C = 1; % Set cluster
    end
    N = 50 * mix.m; % Number of data points
    % Sample data
    data = spm_samp_mix(mix, N);
else
    % If data is a vector, remove NaNs and Inf values
    data = data(isfinite(data));
    mix = 'Own data given';
end

logev = nan(maxcomps, 1);   % Array for log evidence score
vbmix = cell(maxcomps, 1);  % Cell array to store models

for m = 1:maxcomps
    if verbosity > 0
        fprintf('Fitting mixture model with %d components\n', m);
    end

    % Variational Bayesian GMM
    vbmix{m} = spm_mix(data, m, verbosity > 1); 
    if verbosity > 0

        % Plot histogram
        figure
        spm_mix_plot1d (data, vbmix{m})
        title('Fitted model')
    end

    % Store log evidence of the fitted model
    logev(m) = vbmix{m}.fm;
end
logev = logev-min(logev);   % Normalise

if verbosity > 0
    if isstruct(mix)
        % Plot true generative model
        figure;
        spm_mix_plot1d (data, mix);
        title('True generative model')
    end

    % Plot bar chart of log evidence vs number of components
    figure
    bar(logev);
    xlabel('Number of mixture components');
    ylabel('Log Evidence');
end

% If no output arguements requested and verbosity enabled, clear outputs
if nargout == 0 && verbosity > 0 
    clear vbmix
end
