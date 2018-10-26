function addNeighbor(self, id1, id2, c1, c2)
% Add neigbor to neigbors list only if was not added.
%   addNeighbor(self, id1, id2, c1, c2)
% Inputs:
%   self = Frame object.
%   id1, id2 = neigbors labeles.
%   c1, c2 = neigbors COMs.
  if(id1 > id2)
      id_a = id2; id_b = id1;
      c_a = c2; c_b = c1;
  else
      id_a = id1; id_b = id2;
      c_a = c1; c_b = c2;
  end
  % Add the pair to the sorted matrix (sorted by id1 and by id2):
  val = [id_a id_b];
  L = size(self.neighbors, 1);
  [i, alreadyAdded] = Frame.FindDoubleSortedElementPosition(self.neighbors, val);
  if(~alreadyAdded)
      if(i <= 1)
          self.neighbors = [val; self.neighbors];
          self.neighborsDists = [norm(c_a - c_b), self.neighborsDists];
          self.neighborsMidPoints = [mean([c_a; c_b]); self.neighborsMidPoints];
      elseif (i > L)
          self.neighbors = [self.neighbors; val];
          self.neighborsDists = [self.neighborsDists, norm(c_a - c_b)];
          self.neighborsMidPoints = [self.neighborsMidPoints; mean([c_a; c_b])];
      else
        temp1 = self.neighbors; temp2 = self.neighborsDists; temp3 = self.neighborsMidPoints;
        self.neighbors = [temp1(1:i-1,:); val; temp1(i:end,:)];
        self.neighborsDists = [temp2(1:i-1), norm(c_a - c_b), temp2(i:end)];
        self.neighborsMidPoints = [temp3(1:i-1, :); mean([c_a; c_b]); temp3(i:end, :)];
      end
  end
end

