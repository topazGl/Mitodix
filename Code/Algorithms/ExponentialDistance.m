function post = ExponentialDistance(X, gm, c)
%%%
% exponential distance membership, as applied in:
% T. Gilad et al, Fully Unsupervised Symmetry-Based Mitosis Detection in
% Time-Lapse Cell Microscopy, 2018.
%
% post = ExponentialDistance(X, gm, c)
% Inputs:
%   X - input N d-dimentional data points (Nxd).
%   gm - GMM distribution struct:
%       ComponentProportion(Cx1), mu (Cxd), Sigma(dxdxC).
%   c - requested cluster (integer).
%%%
    C = length(gm.ComponentProportion);
    N = size(X, 1);
    dist = zeros(N, C);
    for k=1:C
        P = gm.ComponentProportion(k);
        MU = gm.mu(k,:);
        Sigma = gm.Sigma(:,:,k);
        if(size(Sigma, 1) ~= size(Sigma, 2))
            Sigma = diag(Sigma);
        end

        for nn=1:N
            x = X(nn, :);
            dist(nn,k) = (1/P)*sqrt(det(Sigma))*exp(0.5*(x-MU)*inv(Sigma)*(x-MU)');
        end
    end
    post = zeros(N,1);
    D = 1./dist;
    for nn=1:N
        post(nn) = D(nn,c);
        if(C>1)
            post(nn) = post(nn)/sum(D(nn,:));
        end
    end
    if(C==1)
        post = post/max(post);
    end
end