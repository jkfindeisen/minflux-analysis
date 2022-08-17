function id = minflux_cluster_identify(pos, epsilon, minpts)
% Identifies clusters (by various methods)
%
% Methods
%   - Gaussian rendering and estimating local maxima above a certain
%     threshold
%   - dbscan based (if possible with gaussian mixture afterwards
%   - tesselation based (inspired by Marcel Leutenegger and Henrik von der
%     Emde)

%% parameters default
if nargin < 3 
    epsilon = 10e-9; % dbscan parameter
    minpts = 3;  % dbscan parameter
end

% call to dbscan (this can potentially take very long, scales quite badly
% with number of psotitions)
id = dbscan(pos, epsilon, minpts);

end