function [h, xi, yi, zi, ix, iy, iz] = render_xyz(x, y, z, sx, sy, sz, Rx, Ry, Rz, options)
% Makes a 3D histogram

assert(nargin >= 9);

if nargin < 10 || isempty(options)
    options = struct('type', 'histogram');
end

% get dimensions
Nx = ceil(diff(Rx) / sx);
Ny = ceil(diff(Ry) / sy);
Nz = ceil(diff(Rz) / sz);
h = zeros([Nx, Ny, Nz], 'single');

% get position in pixels
px = (x - Rx(1)) / sx;
py = (y - Ry(1)) / sy;
pz = (z - Rz(1)) / sz;

% get indices
ix = round(px);
iy = round(py);
iz = round(pz);

% remove those outside
m = ix > 0 & ix <= Nx & iy > 0 & iy <= Ny & iz > 0 & iz <= Nz;
px = px(m);
py = py(m);
pz = pz(m);
ix = ix(m);
iy = iy(m);
iz = iz(m);

% switch according to type

switch options.type
    case 'histogram'
        % just a histogram
        
        % fill in histogram (or use accumarray)
        for i = 1 : numel(ix)
            xi = ix(i);
            yi = iy(i);
            zi = iz(i);
            h(xi, yi, zi) = h(xi, yi, zi) + 1;
        end
        
    case 'fixed_gaussian'
        % gaussian with subpixel accuracy
        
        wx = options.fwhm / sx;
        wy = options.fwhm / sy;
        wz = options.fwhm / sz;
        L = ceil(2*max([wx, wy, wz]));
        % small grid
        g = -L:L;
        [xk, yk, zk] = ndgrid(g, g, g);
        
        % remove close at border
        m = ix >= L+1 & ix <= Nx-L-1 & iy >= L+1 & iy <= Ny-L-1 & iz >= L+1 & iz <= Nz-L-1;
        px = px(m);
        py = py(m);
        pz = pz(m);
        ix = ix(m);
        iy = iy(m);
        iz = iz(m);
        
        
        for i = 1 : numel(ix)
            xi = ix(i);
            yi = iy(i);
            zi = iz(i);
            dx = px(i) - xi;
            dy = py(i) - yi;
            dz = pz(i) - zi;
            gx = xi + g;
            gy = yi + g;
            gz = zi + g;
            k = exp(-4*log(2)*((xk-dx).^2/wx^2+(yk-dy).^2/wy^2++(zk-dz).^2/wz^2));
            h(gx, gy, gz) = h(gx, gy, gz) + k;
        end
        
    otherwise
        error('unknown type');
end

% define output xy grid
xi = Rx(1) + (1:Nx) * sx + sx / 2;
yi = Ry(1) + (1:Ny) * sy + sy / 2;
zi = Rz(1) + (1:Nz) * sz + sz / 2;

end