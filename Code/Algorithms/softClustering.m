function [post, gm, c] = softClustering( mitoticFeatures )
% Fuzzy clustering into 2 classes assuming normal distribution. The desired
% class assumed to be dense, with maximal 2nd feature.
% Signature:
%   post = softClustering( mitoticFeatures )
% Inputs:
%   mitoticFeatures: Nx2 features matrix.
% Output:
%   post: posterior probability of the mitotic cluster.
% See also: fitgmdist; posterior
%%
    post = [];
    gm = [];
    c = 1;
    if(isempty(mitoticFeatures))
        return;
    end
    [N, d] =size(mitoticFeatures);
    
    if(N < 5)
        post = ones(N, 1);
    end
    if(d > 2)
        warning('d>2 for 2d clustering');
    end
    d = min(d,2);
    
    mitoticFeatures = mitoticFeatures(:,1:d);
    
    % Expectation Maximization:
    gm = clusterIter(mitoticFeatures);

    % Mitotic Cluster:
    if(d==1)
        [~, c] = min(abs(gm.Sigma));
    else
        % Max symmetry:
        [~, c] = max(gm.mu(:,2));
    end

    % Soft Clustering:
    post = ExponentialDistance(mitoticFeatures, gm, c);
%     post = posterior(gm, mitoticFeatures);
%     post = post(:, c);
end

function [gm] = clusterIter(mitoticFeatures, numIter)
    % Choose cluster with minimal Sigma among the iterations:
    if(~exist('numIter', 'var') || isempty(numIter))
        numIter = 10;
    end
    d = size(mitoticFeatures, 2);
%     inds = (mitoticFeatures(:, 1)>0.1) & (mitoticFeatures(:, 2)>0.5);
    inds = (mitoticFeatures(:, 1)>0.1) & (mitoticFeatures(:, 2)>0.1);
%     inds = (((1-mitoticFeatures(:, 1)) < mitoticFeatures(:, 2)) & (mitoticFeatures(:, 1) > 0) & (mitoticFeatures(:, 1) < 1) & (mitoticFeatures(:, 2) > 0) & (mitoticFeatures(:, 2) < 1));
    N = sum(inds);
    if(N < 100) 
        disp('Not enough candidates for 2 clusters, performing 1 class MLE.');
        C = 1;
        MED = median(mitoticFeatures, 1);
        MAD = median(abs(mitoticFeatures-repmat(MED, [size(mitoticFeatures, 1), 1])));
        mitoticFeatures = mitoticFeatures((abs(mitoticFeatures(:,1)-MED(1)) < 4*MAD(1)) & (abs(mitoticFeatures(:,2)-MED(2))<3*MAD(2)), :);
        if(N < 5)
            if(d==2)
                mu = [median(mitoticFeatures(:,1)), mean(mitoticFeatures(mitoticFeatures(:,2)>mean(mitoticFeatures(:,2)),2))];
                sigma = 0.5*[std(mitoticFeatures(:,1)), 0; 0, std(mitoticFeatures(:,2))];
            else
                mu = [median(mitoticFeatures(:,1)), mean(mitoticFeatures(mitoticFeatures(:,2)>mean(mitoticFeatures(:,2)),2)), median(mitoticFeatures(:,3))];
                sigma = 0.5*[std(mitoticFeatures(:,1)), 0, 0; 0, std(mitoticFeatures(:,2)), 0; 0, 0, std(mitoticFeatures(:,3))];
            end
            gm = gmdistribution(mu,sigma,0.5);
            N = size(mitoticFeatures, 1);
            return;
        end
        numIter = 1;
        N = size(mitoticFeatures, 1);
%         mitoticFeatures = mitoticFeatures(inds, :);
    else
        mitoticFeatures = mitoticFeatures(inds, :);
        C = 2;
    end
    options.MaxIter = 1000;
    
    gm_k = cell(C, 1);
    minSig_k = zeros(C, 1);
    c_K = zeros(C, 1);
    
    for k=1:C
        objs = cell(numIter, 1);
        minSig = zeros(numIter, 1);
    %     symm = zeros(numIter, 1);
    %     symm_sig = zeros(numIter, 1);
        for n=1:numIter
    %         if(n == 1)
    %             lbls_guess = ((mitoticFeatures(:,2) > 0.85) & (abs(mitoticFeatures(:,1) - 0.5) < 0.05)) + 1;
    %             if(C == 1)
    %                 S = [];
    %                 S.Sigma = std(mitoticFeatures(lbls_guess == 2, :));
    %                 S.mu = mean(mitoticFeatures(lbls_guess == 2, :));
    %                 objs{n} = fitgmdist(mitoticFeatures,1, 'Options', options, 'start', S, 'CovarianceType', 'diagonal');
    %             else
    %                 objs{n} = fitgmdist(mitoticFeatures,k, 'Options', options, 'start', lbls_guess);
    %             end
    %         else
                try
                    objs{n} = fitgmdist(mitoticFeatures,k, 'Options', options, 'start', 'plus', 'CovarianceType', 'diagonal');
                catch
                    minSig(n) = inf;
                    continue;
                end
    %         end

            [~,mitotic_cluster] = max(objs{n}.mu(:,2));
            sigma = objs{n}.Sigma(:,:,mitotic_cluster);
%             b = abs(objs{n}.Sigma(:,:,mitotic_cluster));
%             minSig(n) = max(b(:)); 
            if(size(sigma, 1) == size(sigma, 2))
                minSig(n) = abs(sigma(2,2));
            else
                minSig(n) = abs(sigma(2));
            end
            
        end

        [minSig_k(k), c_K(k)] = min(minSig);
        gm_k{k} = objs{c_K(k)};
    end
    [~, k] = min(minSig_k);
    gm = gm_k{k};
end
