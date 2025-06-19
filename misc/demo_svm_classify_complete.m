
clear all
close all

% This script demonstrates how to load a dataset and check the
% classification performance

% Choose dataset
data = 'iono';

switch data
    case 'iono',
        load ionosphere % Load X (features) and Y (labels)
        R = svm_classify(X,Y)   % SVM classifier

    case 'iris',
        load fisheriris % Load 'meas' (features) and 'species' (labels)
        x = meas;
        y = species;
        R = svm_classify(x,y)   % SVM classifier
end
