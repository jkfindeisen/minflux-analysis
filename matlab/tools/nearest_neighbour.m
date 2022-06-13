function [nn] = nearest_neighbour(pos, frames)
% Computes the next neighbour distances (and indices) for each position of
% a (x1, x2, x3, ...) list.
%
% Syntax
%   function [nn] = omex_nearest_neighbour(pos, frames)
%
% Input parameter
%   pos     a (x1, x2, ...,xn) position list (Nxn)
%   frame   [optional] a frame list with equal number of rows as parameter
%           pos indicating that certain objects belong to a certain group
%           and distances are only computed between position with the same
%           frame number (Nx1)
%
% Output parameter
%   nn      a vector with nearest neighbour distances and indices of the
%           particle which is a nearest neighbour (or 0 if there is none)

% error checks
if nargin < 1
    error('Not enough arguments!');
end

if nargin < 2
    frames = ones(size(pos, 1), 1);
end

if size(pos, 1) ~= size(frames, 1)
    error('Parameter pos and frames should have same number of rows!');
end

% initialization

% number of particles
np = size(pos, 1);
nn = zeros(np, 2);
if np == 0
    return;
end
% number of different frames
nf = unique(frames)';

% loop over the frames
for kf = nf
    % get the positions of the actual frame
    idx = find(frames == kf);
    pos2 = pos(idx, :);
    
    % for each object in this frame
    np2 = size(pos2, 1);
    nn2 = zeros(np2, 2);
    for ki = 1 : np2
        % get a single position
        p = repmat(pos2(ki, :), np2, 1);
        
        % compute distance
        d = sqrt(sum((pos2 - p).^2, 2));
        
        % make self-distance solution the least favorable
        d(ki) = max(d) + 1E3;
        
        % find minimum add to output
        [dmin, id] = min(d);
        
        % add to output
        nn2(ki, :) = [dmin, id];
    end
    
    % order output
    nn2(:,2) = idx(nn2(:,2));
    nn(idx, :) = nn2;
end

end