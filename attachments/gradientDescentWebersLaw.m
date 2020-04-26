
clear all;clc;close all;
%% Initialization %%%
dim = 3;
numVec = 10000;
step_size = 0.001;
maxIter = 50;
maxLoss = 10^(-10);
%Ae = generateSPDmatrix(dim);
Z = zeros(dim,numVec);
error = zeros(numVec,1);

%% Input %% 
 X = rand(dim,numVec);
 Y = rand(dim,numVec);
 Z = Y-X;
 %mu = rand(1,dim);
 %mu = [0,0,0];
 %sigma = rand(dim,dim);
 %sigma = eye(3);
 %sigma = sigma*sigma';
 %Z = mvnrnd(mu, sigma, numVec);
 %Z = Z';
Z(1,:) = 2*abs(Z(1,:));
Z(2,:) = abs(Z(2,:));
Z(3,:) = 5*abs(Z(3,:));
%save('X10000.mat', 'X')
%load('X10000.mat');
%save('Y10000.mat', 'Y')
%load('Y10000.mat');
%save('Z10000.mat', 'Z')
%load('Z10000.mat');

error = computedError(X, Z);
%% gradient descent step %%
trainingSize = 9000;
trainingZ = Z(:, 1:trainingSize);
actualTrainingError = error(1:trainingSize);
%[optValGradientDescent,A] = gradientDescentFunc( trainingZ, actualTrainingError, maxIter, maxLoss , step_size);

% Optimal Value %%
[OptTrainingValue, M_closedForm] = closedForm( trainingZ, actualTrainingError);
M_closedForm
OptTrainingValue
%% Testing %%
testingZ = Z(:, trainingSize+1:numVec);
estimatedTestingError = zeros ((numVec - trainingSize),1);
estimatedTestingError = diag(testingZ'*M_closedForm*testingZ);
%estimatedTestingError = tanh(diag(testingZ'*A*testingZ));
diffError = estimatedTestingError - error(trainingSize+1:numVec);
optTestingValue = (1/(numVec - trainingSize))*sum(diffError.^2);

%% Plot error %%
magX = (norm2Signal(X));
magY = (norm2Signal(Y));
magZ = (norm2Signal(Z));
% testingMagX = magX(trainingSize+1:numVec);
% testingMagY = magY(trainingSize+1:numVec);
% testingMagZ = magZ(trainingSize+1:numVec);
% 
actualTestingError = error(trainingSize+1:numVec);
figure(1);
stem(estimatedTestingError,actualTestingError);
xlabel('Estimated d for testing signals');
ylabel('C for testing signals ');
% figure(2);
% 
% surf(Z(1,:), Z(2,:), Z(3,:), estimatedTestingError);
% scatter(testingMagX(actualTestingError == 1), testingMagZ(actualTestingError == 1), 10, 'b','*')
% %scatter(testingMagX(actualTestingError > 0.05), testingMagZ(actualTestingError > 0.05), 10, 'b','*')
% xlabel('|x|');
% ylabel('|z|');
% title('blue points correspond to signal with c = 1');
% hold on;
% scatter(testingMagX(actualTestingError == 0), testingMagZ(actualTestingError == 0), 10, 'r','*')
% %scatter(testingMagX(actualTestingError < 0.05), testingMagZ(actualTestingError < 0.05), 10, 'r','*')
% xlabel('|x|');
% ylabel('|z|');
% figure(3);
% scatter(testingMagX, testingMagZ, 50,estimatedTestingError,'filled');
% xlabel('|x|');
% ylabel('|z|');
% title('Blue points correspond to signal with small estimated d = z^T*M*z');
% figure(4);
% scatter(testingMagX(estimatedTestingError > 0.05), testingMagZ(estimatedTestingError > 0.05), 10, 'b','*');
% xlabel('|x|');
% ylabel('|z|');
% title('Threshold  d = 0.05');
% hold on;
% scatter(testingMagX(estimatedTestingError < 0.05), testingMagZ(estimatedTestingError < 0.05), 10, 'r','*');
% xlabel('|x|');
% ylabel('|z|');
% title('Threshold  d = 0.5');










