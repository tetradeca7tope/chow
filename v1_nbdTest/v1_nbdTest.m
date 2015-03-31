% Estimate the neighborhood of a given node

clear all;
close all;
addpath /usr0/home/kkandasa/projects/Robins/if-estimators/
addpath ../utils/
addpath ../MR

% Problem params
p = 5; d = 3; n = 3000;
% p = 12; d = 4; n = 10000;
nodes = 1:p;

test = 'gauss';
% test = 'mrf';

% Generate Data
switch test
  case 'gauss',  G = GaussGraph(p, d);
  case 'mrf', G = MRFGraph(p, d);
end
X = G.sample(n);

if strcmp(test, 'gauss'), S = X'*X /n; S, inv(G.invSigma), end; % A sanity test

% print test info 
fprintf('Experiment %s: (p,d)=(%d,%d), n = %d\n', test, p, d, n);
G.printNeighbors();

% Prep args for learnNbd
estimCondMIParams.kdePickMethod = 'cv';
estimCondMI = @(S1, S2, S3) condShannonMI(S1, S2, S3, [], []);
params.maxNumNeighbors = 3*d;
params.toPrune = true;
params.thresholdAdd = -inf;
params.thresholdPrune = -inf;

% Call learnNbd
for node = nodes
  neighborsLearned = learnNbd(X, node, estimCondMI, params);

  % Compare neighbors
  fprintf('\nNode: %d, Learned: %s, True: %s\n\n', node, ...
    mat2str(neighborsLearned), mat2str(G.neighbors{node}) );

end

