function phi = CalcPHI(sdf, epsilon)
% Inputs:
% 	sdf: a sign distance function of an image.
%   epsilon: = certainty: 0<epsilon (optional. default is pi).
% Output:
%   phi: sigmoid of the input SDF.
    if(nargin < 2)
        epsilon = pi;
    end
    if(epsilon < 0)
        throw('Input epsilon cannot be negative.');
    end
    if(any(isnan(sdf)))
        [M, N] = size(sdf);
        phi = double(ones(M, N))/(M*N);
    else
        phi = 1./(1+exp(-2*double(sdf)./epsilon));
    end
end

