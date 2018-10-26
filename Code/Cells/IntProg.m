classdef IntProg < handle
    %INTPROG - this class holds data and functions to solve an integer
    %programming problem.
    
    properties(Constant)
        Log = SingleInstance.Logger;        % Application logger.
    end
    
    properties (GetAccess = private)
        daughters_ids = [];
        mothers_ids = [];
        candidates = [];
        Kd = 0;
        Km = 0;
    end
    
    properties (SetAccess = private)
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        M = [];
        f = [];
        N = 0;
    end
    
    properties (Dependent)
        K;
    end
    
    methods
        % Constructor
        function self = IntProg(candidates, NumOfChosen, eqFlag)
        % constructor
        %   self = IntProg(candidates, NumOfChosen, eqFlag)
            try
                if(~exist('eqFlag', 'var') || isempty(eqFlag))
                    eqFlag = 0;
                end

                self.candidates = candidates;
                C = length(candidates);

                % build IDs list:
                daughtersIds = [];
                mothersIds = [];
                for c=1:C
                    cand = candidates(c);
                    daughtersLbls = [cand.Daughter1Label, cand.Daughter2Label];
                    motherLbl = cand.MotherLabel;

                    daughtersIds = [daughtersIds, daughtersLbls];
                    mothersIds = [mothersIds, motherLbl];
                end
                self.daughters_ids = sort(unique(daughtersIds));
                self.mothers_ids = sort(unique(mothersIds));

                self.Kd = length(self.daughters_ids);
                self.Km = length(self.mothers_ids);
                self.N = C;%sum(1:self.K);
                self.b = ones(self.K, 1);
                self.f = 10*ones(self.N, 1);

                % build f and A:
                self.A = zeros(self.K, self.N);
                % Each col in A is a candidate. The rows are all the duaghters
                % followed by all the mothers. Each col is zeros but to 3
                % values: 2 daughters and the mother are set to '1'.
                % f is a vector of scores for each candidate.
                r = 1;
                for c = 1:C
                    pair = [candidates(c).Daughter1Label,  candidates(c).Daughter2Label];
                    mother = candidates(c).MotherLabel;
                    score = candidates(c).MitosisScore;

                    i = find(self.daughters_ids==pair(1), 1, 'first');
                    j = find(self.daughters_ids==pair(2), 1, 'first');
                    k = self.Kd + find(self.mothers_ids==mother, 1, 'first');
                    self.A(i,r) = 1; self.A(j,r) = 1; self.A(k,r) = 1;
                    self.f(r) = -1*score;
                    r = r+1;
                end

                % constraint number of chosen elements:
                if(exist('NumOfChosen', 'var') && ~isempty(NumOfChosen))
                    if(eqFlag)
                        self.Aeq = ones(1, C);
                        self.beq = NumOfChosen;
                    else
                        self.A = [self.A; ones(1, C)];
                        self.b = [self.b; NumOfChosen];
                    end
                end
            catch EX
                rethrow EX
            end
        end
        
        % Get K
        function K = get.K(self)
            K = self.Kd + self.Km;
        end
        
        %% Public methods
        
        % solve integer programming problem
        function [chosen, idxs, scores] = solve(self)
        % [chosen, idxs, scores] = solve(self)
            options = optimoptions('intlinprog','Display','off');
            x = intlinprog(self.f, 1:self.N, self.A, self.b, self.Aeq, self.beq, zeros(1, self.N), ones(1, self.N), options);
            idxs = find(x~=0);
            if (isempty(idxs))
                chosen = [];
                return;
            end
                
            chosen = self.candidates(idxs);
            scores = (-1)*self.f(idxs);
%             L = length(idxs);
%             pairs = zeros(L,2);
%             for i=1:L
%                 n = idxs(i);
%                 p = find(self.A(n,:)~=0);
%                 if(length(p) == 2)
%                     pairs(n,:) = [self.ids(p(1)), self.ids(p(2))];
%                 end
%             end
        end
    end
    
end

