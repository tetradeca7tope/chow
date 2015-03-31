function utLearnNbd
% Unit Test for learnNbd

  % Test 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('Test 1: Mock Example\n'); 
  fprintf('================================================================\n');
  X = repmat(1:20, 500, 1);
  estimCondMI = @mockEstimCondMI;
  params.thresholdAdd = 1;
%     % Test with different parameters
%     params.maxNumNeighbors = 20; params.toPrune = true;
%     params.thresholdPrune = 1.15; neigbors = learnNbd(X, 1, params);
    % Test with different parameters
    params.maxNumNeighbors = 5; params.toPrune = false;
    params.thresholdPrune = 1.15; neigbors = learnNbd(X, 1, estimCondMI,params);

  fprintf('Super set of True Neighbors: 2, 5, 8, 11, 14, 17\n');
  fprintf('Learned: %s\n', mat2str(neigbors));
  

end


function val = mockEstimCondMI(X, Y, Z)

  x = X(1, :);
  y = Y(1, :);
  if ismember(x+y, 3:3:20)
    val = 1 + 0.5*rand();
  else
    val = 0 + 0.5*rand();
  end

end

