function [pos, id, rpos] = minflux_filter(pos, id, epsilon, minpts)
% Filters localizations
%
% Note
% - Filtering by minimal number of points in dbscan follows an idea of
%   Mariano Bossi. 
%
% Other ideas include dbscan, density based.

%% parameters default
if nargin < 3 
    epsilon = 10e-9; % dbscan parameter THIS IS 10 nm!
    minpts = 3;  % dbscan parameter
end

%% sort by id and obtain first, last indices (id should already be sorted, but you never know)
[id, idx] = sort(id);
pos = pos(idx, :);
i1 = find(diff([0; id]) > 0);
i2 = find(diff([id; Inf]) > 0);

%% loop over all unique ids, call dbscan and create new ids
k = 1;
id(:) = 0;
rpos = zeros(size(pos));
for i = 1 : numel(i1)
    g = i1(i) : i2(i);
    pos_i = pos(g, :);
    
    % call to dbscan
    idx = dbscan(pos_i, epsilon, minpts);    
    
    % all that are inside a cluster, write out with new cluster ids
    uidx = unique(idx(idx > 0));
    for j = 1 : numel(uidx)
        ix = idx == uidx(j);
        gj = g(ix);
        % add new cluster index
        id(gj) = k;
        k = k + 1;
        % get rpos = pos - center of cluster
        pj = pos_i(ix, :);
        rpos(gj, :) = pj - mean(pj, 1);
    end
end

% keep only those within clusters
idx = id > 0;
pos = pos(idx, :);
id = id(idx);
rpos = rpos(idx, :);

end