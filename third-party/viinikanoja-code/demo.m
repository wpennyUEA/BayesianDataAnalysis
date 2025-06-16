%
% Simple illustration of the robustness effect on the conditional means
%
% Copyright (C) 2010 Jaakko Viinikanoja
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% Depends on Lightspeed toolbox by Tom Minka
%

%addpath('./lightspeed/');


N = 100;
M = 1;

dX = 2;
dY = 2;
D = 2;

latentDim = 2;

clear nModel;
clear tModel;

% Generate some toydata
dataX = randnorm(N, 0.2*ones(dX,1), [], eye(dX));

% ypart
W1 = [1 0; 1 0] + 0.1*rand(2);

% Noise matrix cholesky factor
noiseY = diag(0.2*ones(2,1));

%dataY(1,:) = 
dataY = W1*dataX + randnorm(N, zeros(2,1),[], noiseY);

% Outliers from uniform distributions
NOutliers = 5;
outliersX = (rand(2,NOutliers)-0.5)*10;
outliersY = (rand(2,NOutliers)-0.5)*10;

%Datasets with outliers
oDataX = [dataX outliersX];
oDataY = [dataY outliersY];

% Learning
%---------------------------------------------
learningParameters = vbmcca.initLearningParameters(150, 10^(-6));

% Initialise model with normal latents
nModel = vbmcca(M,dX,dY,D);
nModel.normalLatents = 1;

% Robust model
tModel = vbmcca(M,dX,dY,D);

nModel = nModel.initWithKMeans(oDataX, oDataY);
tModel = tModel.initWithKMeans(oDataX, oDataY);

[nModel, nEnergies] = nModel.learn(oDataX, oDataY, learningParameters); 
[tModel, tEnergies] = tModel.learn(oDataX, oDataY, learningParameters); 

% Predict new unobserved values
%-------------------------------------------------
newN = 8;
newdataX = randnorm(newN, zeros(dX,1), [], eye(dX));
newdataY = W1*newdataX;

[eYn,covYn] = nModel.inferConditionalYMean(newdataX);
[eYt,covYt] = tModel.inferConditionalYMean(newdataX);

%Prediction error
%---------------------------------------
nError = 0;
tError = 0;
for i = 1 : newN
    nError = nError + sqrt( sum((eYn(:,i) - newdataY(:,i)).*(eYn(:,i) - newdataY(:,i))) );
    tError = tError + sqrt( sum((eYt(:,i) - newdataY(:,i)).*(eYt(:,i) - newdataY(:,i))) );
end

disp(strcat('Total prediction error(n): ',num2str(nError)));
disp(strcat('Total prediction error(t): ',num2str(tError)));