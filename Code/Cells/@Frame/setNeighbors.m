function setNeighbors( self )
% Calculate neigbors graph in frame.
% Create a graph using Delaunay Triangulation. Each cell center is a node.
%   [ neighbors, dists, midPoints ] = calcNeighbors( self, centers, labels )
% Inputs:
%   self = Frame object.

    self.neighbors = [];       % An 2xM matrix of the IDs of all connected nodes.
    self.neighborsDists = [];  % An 1xM vector consits of the length of each edge
                               % (Euclidian distance between the center points of the two
                               % connected nodes).
    self.neighborsMidPoints = []; % An 2xM matrix of the mid points of connected nodes.
    
    % Check that input was not empty:
    if(isempty(self.labels))
        return;
    end

    % Get neighbors according to Delaunay Triangolation graph
    % connections:
    try
        L = min(max(self.CellLength.STD, 0.25*self.CellLength.Mean), 1);
        R = 0.5*min([size(self.Image, 1), size(self.Image, 2)]);%5*self.CellLength.Mean;
        [DT_connect, connections, connect, iterCount, converg] = StochasticDelaunay(self.centers, L, 'normal', R);%'uniform', R);
        
    catch EX
        self.Log.mlog = {LogType.ERROR ...
                      ,mfilename('class') ...
                      ,[self.Log.Me,' Problem using delaunayTriangulation.']};
        rethrow(EX);
    end
% 
%     im = self.Image;
%     mx = min(im(:));
%     Mx = max(im(:));
%     im = im - mx;
%     im = im / Mx;
%     f = figure; h = imagesc(im); hold on; ax = ancestor(h, 'axes'); colormap gray;
%     triplot(DT_connect, self.centers(:,1), self.centers(:,2), 'r', 'Parent', ax);
%     triplot(connect, self.centers(:,1), self.centers(:,2), 'c', 'Parent', ax);
%     hold off;
%     close(f);
%     disp(iterCount);
if(iterCount>15)
    f = figure; plot(converg, '*'); title(['#Iterations = ', num2str(iterCount)]);
    close(f);
end

    N = size(connections,1);
%     N = floor(N/2); % matrix is symmetric
    for i=1:N
          % For each triangle, create 3 pairs of neighbors:
          idx1 = i;
          id1 = self.labels(idx1);
          c1 = self.centers(idx1,:);
          
          idx = find(connections(idx1,:) == 1);
          idx = idx(idx > i); % Matrix is symmetric
          M = length(idx);
          for m=1:M
              idx2 = idx(m);
              c2 = self.centers(idx2,:);
              id2 = self.labels(idx2);
              self.addNeighbor(id1, id2, c1, c2);
          end
    end
end

