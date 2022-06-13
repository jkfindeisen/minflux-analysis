function idx = density_filter(pos, ROI, fwhm, T)
% Density based filtering
%   pos     Nx2 or Nx3 positions
%   ROI     2xD ranges
%   fwhm    FWHM for blurring
%   T       threshold after blurring
%
%   idx     true if should be kept, false if should be discarded

assert(nargin == 4);

px = min(fwhm) / 3;

D = size(pos, 2); % dimensionality

% rendering with Gaussian FWHM
switch D
    case 2
        [h, ~, ~, idx, m] = render_xy(pos(:, 1), pos(:, 2), px, px, ROI(:, 1), ROI(:, 2), struct('type', 'fixed_gaussian', 'fwhm', fwhm));
    case 3
        [h, ~, ~, ix, iy, iz, m] = render_xyz(pos(:, 1), pos(:, 2), pos(:, 3), px, px, px, ROI(:, 1), ROI(:, 2), ROI(:, 3), struct('type', 'fixed_gaussian', 'fwhm', fwhm));        
        idx = sub2ind(size(h), ix, iy, iz);
    otherwise
        error('Unsupported dimensionality %d', D);
end

% density based filtering, those outside of ROI will also be filtered away
keep = h > T;
keep = keep(idx);

idx = false([size(pos, 1), 1]);
ix = find(m);
idx(ix(keep)) = true;

end