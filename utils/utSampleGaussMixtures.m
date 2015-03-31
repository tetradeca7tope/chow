% Unit test for sampleGaussMixtures

% close all;
clear all;

numCentresCandidates = 2:5;
numSamples = 1000;
numDims = 1;

for numCentres = numCentresCandidates
  subplot(131);
  X = sampleGaussMixtures(numSamples, 1, numCentres);
  plot(X, rand(numSamples, 1), 'kx');

  subplot(132);
  X = sampleGaussMixtures(numSamples, 2, numCentres);
  plot(X(:, 1), X(:,2), 'kx');

  subplot(133);
  X = sampleGaussMixtures(numSamples, 3, numCentres);
  plot3(X(:, 1), X(:,2), X(:,3), 'kx');

  pause,
end

