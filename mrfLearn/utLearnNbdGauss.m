% Tests the Algorithm with the true mutual Informations of Gaussians
function utLearnNbdGauss 

  clear all;
  close all;

  % Problme parameters
  p = 5; d = 3; 
  p = 40; d = 7;
  n = 100;
  nodes = 1:p;

  % Generate the Graph
  G = GaussGraph(p, d);
  data = repmat(1:p, n, 1);
  C = inv(G.invSigma);

  % Prepare arguments for learnnbd
  estimCondMI = @(X, Y, Z) mockEstimCondMI(X, Y, Z, C); 
  params.maxNumNeighbors = 3*d;
  params.toPrune = true;
  params.thresholdAdd = -inf;
  params.thresholdPrune = -inf;

  % Call learnNbd
  for node = nodes
    neighborsLearned = learnNbd(data, node, estimCondMI, params);
          
    % Compare neighbors
    fprintf('\nNode: %d, Learned: %s, True: %s\n\n', node, ...
      mat2str(neighborsLearned), mat2str(G.neighbors{node}) );
          
  end

end


function val = mockEstimCondMI(X, Y, Z, C)

  x = X(1, :);
  y = Y(1, :);
  z = Z(1, :);

  Cxyz = C([x, y, z], [x, y, z]);
  Cz = C(z, z);
  Cxz = C([x, z], [x, z]);
  Cyz = C([y, z], [y, z]);

  val = -0.5*( log(det(Cxyz)) + log(det(Cz)) - log(det(Cxz)) - log(det(Cyz)) );

end

