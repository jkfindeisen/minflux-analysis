function [h, xi, yi, idx, m] = render_xy(x, y, sx, sy, Rx, Ry, options)
% Makes a 2D histogram. Just to stay compatible with render_xyz

assert(nargin >= 6);

if nargin < 7 || isempty(options)
    options = struct('type', 'histogram');
end

% get dimensions and allocate output array
Nx = ceil(diff(Rx) / sx);
Ny = ceil(diff(Ry) / sy);
h = zeros([Nx, Ny], 'single');

% get position in pixels
px = (x - Rx(1)) / sx;
py = (y - Ry(1)) / sy;

% get position indices
ix = round(px);
iy = round(py);

% remove those outside
m = ix > 0 & ix <= Nx & iy > 0 & iy <= Ny;
px = px(m);
py = py(m);
ix = ix(m);
iy = iy(m);

% get linear output index
idx = ix + (iy-1)*Nx;

% switch according to type

switch options.type
    case 'histogram'
        % just a histogram
        
        % fill in histogram (or use accumarray)
        for i = 1 : numel(ix)
            xi = ix(i);
            yi = iy(i);
            h(xi, yi) = h(xi, yi) + 1;
        end
        
    case 'fixed_gaussian'
        % gaussian with subpixel accuracy
        
        wx = options.fwhm / sx;
        wy = options.fwhm / sy;
        L = ceil(2*max(wx, wy));
        % small grid
        g = -L:L;
        [xk, yk] = ndgrid(g, g);
        
        % remove close at border
        m = ix >= L+1 & ix <= Nx-L-1 & iy >= L+1 & iy <= Ny-L-1;
        px = px(m);
        py = py(m);
        ix = ix(m);
        iy = iy(m);
        
        for i = 1 : numel(ix)
            xi = ix(i);
            yi = iy(i);
            dx = px(i) - xi;
            dy = py(i) - yi;
            gx = xi + g;
            gy = yi + g;
            k = exp(-4*log(2)*((xk-dx).^2/wx^2+(yk-dy).^2/wy^2));
            h(gx, gy) = h(gx, gy) + k;
%             h(gx, gy) = max(h(gx, gy), k); % this would 
        end
        
    otherwise
        error('unknown type');
end

% define output xy grid
xi = Rx(1) + (1:Nx) * sx + sx / 2;
yi = Ry(1) + (1:Ny) * sy + sy / 2;

end