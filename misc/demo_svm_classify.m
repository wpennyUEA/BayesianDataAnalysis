
clear all
close all

data = 'iono';

switch data
    case 'iono',
        load ionosphere
        R = svm_classify(X,Y)

    case 'iris',
        load fisheriris
        x = meas;
        y = species;
        R = svm_classify(x,y)
end
