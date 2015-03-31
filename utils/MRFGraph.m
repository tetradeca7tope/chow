classdef MRFGraph < handle

  properties
    p; % number of nodes
    d; % maximum degree
    nodes; % The nodes in the graph
    neighbors; % The neighbors of each node
    parents; % The parents of each node
    children; % The children of each node
    atoms; % The atoms - one per parent. To be used for sampling.
    coeffs; % The coeffs - one per parent. To be used for sampling.
    numBaseNodes; % Number of base nodes.
    numChildren; % A p-vector to keep track of the number of children
    numParents; % A p-vector to keep track of the number of parents
    numNeighbors; % A p-vector to keep track of the number of neighbors
    atomCandidates = {@(t)t,  @(t)2.1*t,  @(t)sign(t).*abs(t).^1.1,  ...
      @(t)sign(t).*abs(t).^0.9, @(t)sign(t).*abs(t).^0.8};
    addNoiseStrength = 0.001;
    numBaseGaussMixtures = 2; % Number of base mixtures to sample from
  end

  methods

    % Constructor
    function obj = MRFGraph(p, d)
      obj.p = p;
      obj.d = d;
      obj.nodes = [1:p]';
      obj.neighbors = cell(p,1);
      obj.parents = cell(p,1);
      obj.children = cell(p,1);
      obj.numBaseNodes = max(d+1, min(floor(p/2), 2*d)); 
      obj.numNeighbors = zeros(p, 1);
      % Now generate the graph
      obj.genGraph();
    end

    % Function which generates graph
    function genGraph(obj)

      % Prelims
      d = obj.d;
      p = obj.p;
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
        % Determine the coefficients
        parPreWeights = randn(obj.numParents(i), 1);
        obj.coeffs{i} = parPreWeights / sum(abs(parPreWeights));
        obj.atoms{i} = cell(obj.numParents(i), 1);

        % Now process each parent
        for j = 1:obj.numParents(i)
          currParent = obj.parents{i}(j);
          parCandIdx = parCandIdxs(j);
          
          % Determine the atom
          atomIdx = randi([1,numel(obj.atomCandidates)]);
          obj.atoms{i}{j} = obj.atomCandidates{atomIdx};

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

    end

    % Function to sample data from this graph
    function X = sample(obj, numSamples)
      % First generate Gaussians
      X = sampleGaussMixtures(numSamples, obj.p, obj.numBaseGaussMixtures);
      % Now go through each node in order. The nodes are already in topological
      % order so no need to keep checking back.
      for i = 1:obj.p
        if ~isempty(obj.parents{i})
          Xcurr = zeros(numSamples, 1);
          for j = 1:numel(obj.parents{i})
            Xcurr = Xcurr + ...
                    obj.coeffs{i}(j) * obj.atoms{i}{j}( X(:,obj.parents{i}(j)));
          end
          X(:,i) = Xcurr + obj.addNoiseStrength*rand(numSamples, 1);
        end
      end
    end

    function printChildren(obj)
      obj.printGraphAttribute(obj.children, 'Children');
    end

    function printParents(obj)
      obj.printGraphAttribute(obj.parents, 'Parents');
    end

    function printNeighbors(obj)
      obj.printGraphAttribute(obj.neighbors, 'Neighbors');
    end

    function printGraphAttribute(obj, attr, attrStr)
      fprintf('%s\n', attrStr);
      for i = 1:obj.p
        fprintf('%d: (%d) %s\n',i, numel(attr{i}), mat2str(attr{i}));
      end
      fprintf('\n');
    end

    % Prints out neighbors, parents and children.
    function printGraph(obj)
      fprintf('# Nodes: %d, Max-deg: %d, NumBaseNodes: %d,  ', ...
        obj.p, obj.d);
      obj.printParents(); 
      obj.printChildren(); 
      obj.printNeighbors(); 
    end

    % Prints out the transformations of each child as a function of its parent.
    function printMRFStructure(obj)
      fprintf('# Nodes: %d, Max-deg: %d, NumBaseNodes: %d,  ', ...
        obj.p, obj.d);
      for i = 1:obj.p
        printStr = sprintf('%d: ', i);
        for j = 1:obj.numParents(i)
          printStr = sprintf('%s (%d, %0.4f, %s), ', printStr, ...
            obj.parents{i}(j), obj.coeffs{i}(j), func2str(obj.atoms{i}{j}) );
        end
        fprintf('%s\n', printStr);
      end
      fprintf('\n');
    end

    % Tests if the sample is from this MRF.
    function isSampleMember = testSampleMembership(obj, X)
      % Prelims
      TOL = 1e-10;
      isSampleMember = true;
      n = size(X, 1);

      for i = (obj.numBaseNodes+1):obj.p
        XCurr = X(:, i); 
        XRegen = zeros(n, 1);
        for j = 1:obj.numParents(i)
          XRegen = XRegen + ...
            obj.coeffs{i}(j) * obj.atoms{i}{j}(X(:, obj.parents{i}(j)));
        end
        if ~all(XCurr - XRegen < obj.addNoiseStrength + TOL)
          isSampleMember = false;
          break;
        end
      end
    end

  end % end methods

end % end classdef

