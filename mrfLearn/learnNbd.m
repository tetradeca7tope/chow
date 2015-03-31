function neighbors = learnNbd(X, node, estimCondMI, params)
% Learns the neighborhood of the node given in node (an index \in [1, 2, ... p])
% where p = size(X, 2) is the dimension of the problem.
% X is the data matrix and is of size n x p.
% params is a struct and contains some of the following fields,
%   - maxNumNeighbors: The maximum Number of neighbors to add
%   - toPrune: Should there be a pruning phase after adding ?
%   - thresholdAdd: If the Conditional MI is lower than this, stop adding.
%   - thresholdPrune: If the Conditional MI is lower than this, stop adding.
% Returns: neighbors - a vector of indices indicating the neighbors

  VERBOSE = true;

  % Prelims
  [n, p] = size(X);
  neighbors = zeros(1,0); % Create empty matrix
  nbdCount = 0;

  % Adding Phase
  % ============================================================================
  for i = 1:params.maxNumNeighbors
    
    nonNeighbors = setdiff(1:p, [neighbors; node]);
    numNonNeighbors = numel(nonNeighbors);

    % Obtain the conditional MI with every other non-neighbor
    nonNbdCondMIs = zeros(numNonNeighbors, 1);
    for j = 1:numNonNeighbors
      currNonNb = nonNeighbors(j);
      nonNbdCondMIs(j) = estimCondMI(X(:,node), X(:,currNonNb), X(:,neighbors));
    end
    % Pick the maximum
    [maxCMIVal, maxCMIIdx] = max(nonNbdCondMIs);
        if VERBOSE
          fprintf('Nb#%d: %d, %0.5f\n', i, nonNeighbors(maxCMIIdx), maxCMIVal);
        end

    if maxCMIVal > params.thresholdAdd
      neighbors = [neighbors; nonNeighbors(maxCMIIdx)];
    else
      break;
    end
  end

  % Pruning Phase
  % ============================================================================
  if params.toPrune
    neighborsPruned = neighbors;
          if VERBOSE, fprintf('\n'); end
    for j = 1:numel(neighbors)
      currNb = neighbors(j);
      remNbd = setdiff(neighborsPruned, [currNb]);
      currCMI = estimCondMI(X(:,node), X(:,currNb), X(:,remNbd));
          if VERBOSE
            fprintf('Nb#%d: %d, %0.5f\n', j, neighbors(j), currCMI);
          end
      if currCMI < params.thresholdPrune
        neighborsPruned = remNbd;
      end
    end

    neighbors = neighborsPruned;

  end

end

