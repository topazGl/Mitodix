function [connect, connections, connect_org, iterCount, converg] = StochasticDelaunay(centers, L, noiseType, R, stabilityThr, maxIter)
% Create stochastic Delaunay graph (additional edges for ambigues cases).
%   [connect, connections, connect_org, iterCount, converg] = StochasticDelaunay(centers, L, noiseType, R, stabilityThr, maxIter)
% Input:
%   centers = Nxd graph nodes (N = number of points, d = dimantion).
%   L = maximal possible noise to add to centers (Default=0, output in this case is regular Delaunay).
%   noiseType =  'uniform' for U(0, L) noise or 'normal' for abs of white
%               Gaussian noise truncated between 0 to L  (Default = 'uniform').
%   R = maximal distance between connected edges (Default = inf).
%   stabilityThr = maximal consecutive iterations in which no additional
%                  edges were added to the graph (Default = 3).
%   maxIter = maximal number of iterations (Default = 50).
%
% Output:
%   connect = Mx3 stochastic connect matrix.
%   connections = NxN binary matix, '1' where the nodes are connected in
%                 the output graph.
%   connect_org = Mx3 original Delaunay connect matrix.
%   iterCount = Number of iterations
%   converg = Number of new connections in each iteration.
% Note: this code produces triangles inside other triangles.
%
% Example:
%   rng(100);
%   centers = gallery('uniformdata',40,2,2);
%   [connect, connections, connect_org, iterCount, converg] = StochasticDelaunay(centers, 0.025);
%   figure;
%   triplot(connect, centers(:,1), centers(:,2), 'r'); title('Stochastic Delaunay (new edges in red)');
%   hold on; triplot(connect_org, centers(:,1), centers(:,2), 'c'); hold off;
%
% See also: delaunayTriangulation
% Author: T. Gilad, 2017
    if(~exist('L', 'var') || isempty(L))
                L = 0;
    end
    
    if(~exist('noiseType', 'var') || isempty(noiseType))
            noiseType = 'uniform';
    end
    
    if(~exist('R', 'var') || isempty(R))
                R = inf;
    end
    
    if(~exist('stabilityThr', 'var') || isempty(stabilityThr))
            stabilityThr = 3;
    end
    if(~exist('maxIter', 'var') || isempty(maxIter))
            maxIter = 50;
    end
    
    [connections, ~, connect_org] = DelaunayConnectivityMatrix( centers, [], [], R );
    if(~L)
        return;
    end
    
    stabilityCounter = 0;
    converg = zeros(1, maxIter);
    for j=1:maxIter
        connectionsNoisy = DelaunayConnectivityMatrix( centers, L, noiseType, R );
        newConnections = sum((connectionsNoisy(:) == 1) & (connectionsNoisy(:) > connections(:)));
        converg(j) = newConnections;
        if(newConnections == 0)
            stabilityCounter = stabilityCounter + 1;
        else
            connections = max(connections, connectionsNoisy);
            stabilityCounter = 0;
        end
        if(stabilityCounter == stabilityThr)
            break;
        end
    end
    connect = MatrixToTriangulation(connections);
    iterCount = j;
end

function [connections, centers, connect] = DelaunayConnectivityMatrix( nodes, L, noiseType, R )
% Calculate connectivity matrix according to Delaunay Triangulation (DT).
% Add noise of maximal size L to each of the input nodes and calculate DT
% of the noised nodes. 
%   connections = DelaunayConnectivityMatrix( nodes, L, noiseType, R )
% Inputs:
%   nodes = Nxd graph nodes (N = number of points, d = dimantion).
%   L = (optional) maximal noise size. Noise is truncated between 0 to L.
%        If not given no noise is added.
%   noiseType = (optional) 'uniform' for U(0, L) noise (default) or
%               'normal' for abs of white Gaussian noise with sigma of L/2
%               truncated between 0 to L.
%   R = maximal distance between connected edges (Default = inf).
% Output:
%   connections = NxN binary matix, '1' where the nodes are connected in
%                 the noisy DT.
%   centers = The noisy nodes.
%   connect = Mx3 DT connections matrix (triangles as nodes IDs).
% Author: T. Gilad, 2017
% See also: TruncatedGaussian - https://www.mathworks.com/matlabcentral/fileexchange/23832-truncated-gaussian

    N = size(nodes, 1);
    connections = zeros(N, N);
    
    if(~exist('L', 'var') || isempty(L))
            L = 0;
    end
    
    if(~exist('noiseType', 'var') || isempty(noiseType))
            noiseType = 'uniform';
    end
    noiseType = lower(noiseType);
    
    if(~exist('R', 'var') || isempty(R))
        R = inf;
    end
    
    try
    % Calculate Delaunay Triangulation:
        status = 0;
        counter = 0;
        while(~status && counter < 10)
            centers = addNoise(nodes, L, noiseType);
            DT = delaunayTriangulation(centers);
            if(~isequal(centers, DT.Points))
                counter = counter + 1;
                continue;
            end
            connect = DT.ConnectivityList;
            status = 1;
        end
        if(~status)
            error('Problem in Delaunay Triangulation nodes locations order');
        end
    catch EX
        rethrow(EX);
    end
    
    % Build connectivity matrix:
    M = size(connect,1);
    for i=1:M
          % For each triangle, create 3 pairs of neighbors:
          val = 1;
          if(R<inf)
              d = norm(nodes(connect(i,1),:)-nodes(connect(i,2),:));
              if(d>R)
                  val = 0.5;
              end
          end
          connections(connect(i,1),connect(i,2)) = val;
          connections(connect(i,2),connect(i,1)) = val;
          
          val = 1;
          if(R<inf)
              d = norm(nodes(connect(i,1),:)-nodes(connect(i,3),:));
              if(d>R)
                val = 0.5;
              end
          end
          connections(connect(i,1),connect(i,3)) = val;
          connections(connect(i,3),connect(i,1)) = val;
          
          val = 1;
          if(R<inf)
              d = norm(nodes(connect(i,2),:)-nodes(connect(i,3),:));
              if(d>R)
                  val = 0.5;
              end
          end
          connections(connect(i,2),connect(i,3)) = val;
          connections(connect(i,3),connect(i,2)) = val;
    end
end

function connect = MatrixToTriangulation(connections)
% Convert NxN connectivity matrix to Mx3 Triangulation connect matrix.
%   connect = MatrixToTriangulation(connections)
% Input:
%   connections = NxN binary connections matix.
% Output:
%   connect = Mx3 connect matrix.
% Note: this code produces triangles inside other triangles!
% Author: T. Gilad, 2017

    connect = [];
    N = size(connections, 1);
    for i=1:N
        row = connections(i,:);
        js = find(row == 1);
        for i1=1:length(js)
            j = js(i1);
            col = connections(:,j);
            ks = find(row'.*col == 1);
            for id2 = 1:length(ks)
                k = ks(id2);
                if(connections(i,k) == 1)
                    connect = [connect; [i,j, k]];
                end
            end
        end
    end
    connect = sort(connect, 2);
    connect = unique(connect, 'rows');
end

function noisy_nodes = addNoise(nodes, L, noiseType)
% Generate noise:
	[N, dim] = size(nodes);
    if(isequal(noiseType, 'normal'))
        noise = abs(TruncatedGaussian(L/2, [-L, L], [N, dim]));
    else
        noise = L*rand(N, dim);
    end
    
    % Add noise to nodes:
    noisy_nodes = nodes + noise;
end