classdef GaussGraph < MRFGraph

  properties
    invSigma;
    sigmaChol;
    LT; % generate data via obj.LT\z z ~ N(0, I)
  end

  methods %(Access = public)

    % Constructor
    function obj = GaussGraph(p, d)
      obj = obj@MRFGraph(p, d);
    end

    function genGraph(obj)

      % Prelims
      d = obj.d;
      p = obj.p;

      % We will use Q to generate the inverse covariance matrix
      Q = zeros(p); 

      % Maintain an array of candidate parents and a count of their children
      parentCandidates = (1:obj.numBaseNodes)';
      numChildrenLeft = d*ones(obj.numBaseNodes, 1);

      for i = (obj.numBaseNodes+1):p
        % First determine the parents
        numCandParents = numel(parentCandidates);
        obj.numParents(i) = min(randi(ceil([d/2 d])), numCandParents-1);
        parCandIdxs = sort( datasample(1:numCandParents, obj.numParents(i), ...
          'Replace', false) );
        obj.parents{i} = parentCandidates(parCandIdxs);
        obj.neighbors{i} = obj.parents{i};

        % Now process each parent
        for j = 1:obj.numParents(i)
          currParent = obj.parents{i}(j);
          parCandIdx = parCandIdxs(j);

          % Determine the coefficients in Q
          Q(i, currParent) = (0.25*p + 0.2*p*rand()) * (randi([0 1])*2 -1);
          Q(currParent, i) = Q(i, currParent);

          % Book Keeping
          numChildrenLeft(parCandIdx) = numChildrenLeft(parCandIdx) - 1;
          obj.children{currParent} = [obj.children{currParent}; i];
          obj.neighbors{currParent} = [obj.neighbors{currParent}; i];
        end

        % Add the current node
        parentCandidates = [parentCandidates; i];
        numChildrenLeft = [numChildrenLeft; (d-obj.numParents(i))];
        % Now update parentCandidates
        removeParents = (numChildrenLeft == 0);
        parentCandidates = parentCandidates(~removeParents);
        numChildrenLeft = numChildrenLeft(~removeParents);
      end

      % Now add strong diagonals
      Q(eye(p) == 1) = 1.4*max(max(abs(Q))) + 1;
      [sigmaChol, diagPower] = stableCholesky(Q);
      Q(eye(p) == 1) = Q(eye(p)==1) + 10^diagPower;
      obj.invSigma = Q;
      obj.sigmaChol = sigmaChol;
%       obj.sigmaChol = chol(Q, 'lower');
      obj.LT = obj.sigmaChol'; % generate data via obj.LT\z z ~ N(0, I)

    end

    % Function to override sampling in MRF Graph
    function X = sample(obj, numSamples)
      X = randn(numSamples, obj.p) / obj.sigmaChol;
    end

  end % end methods

end % end classdef
