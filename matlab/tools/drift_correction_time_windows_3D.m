function [dx, dy, dz, dxt, dyt, dzt, ti, fig] = drift_correction_time_windows_3D(x, y, z, t, Rx, Ry, Rz, T, sxyz, use_gpu)
% drift correction, based on auto-correlation
%
% f     frames for each event
% W     chunk size in frames
% x, y, z   coordinates
% Rx, Ry, Rz    ranges of the coordinates
% dxyz      spacing in the coordinate range
% use_gpu   true -> tries to speed up using the gpu
% show_results true -> shows a graph
%
% dx, dy, dz  drift in units of x, y, z

% parameters
if use_gpu
    gpuDevice(1); % resets the memory of the GPU
end

Rt = t([1,end]);
Ns = floor(diff(Rt) / T); % total number of time windows
assert(Ns > 1, 'Need at least two time windows, reduce T');

% parameters
CR = 8; % center of mass
D = Ns; % how many frames maximally in the cross-correlation
w = linspace(1, 0.2, D).'; % weight with distance
% w = ones(D, 1);
l = 0.01;

% get dimensions
Nx = ceil(diff(Rx) / sxyz);
Ny = ceil(diff(Ry) / sxyz);
Nz = ceil(diff(Rz) / sxyz);
c = round([Nx, Ny, Nz]/2) + 1;

% create all the histograms
h = cell(Ns, 1);
ti = zeros(Ns, 1); % average times of the snapshots
opt = struct('type', 'fixed_gaussian', 'fwhm', 3*sxyz);
for j = 1 : Ns
    t0 = Rt(1) + (j-1)*T;
    idx = t >= t0 & t < t0 + T;
    ti(j) = mean(t(idx));
    hj = render_xyz(x(idx), y(idx), z(idx), sxyz, sxyz, sxyz, Rx, Ry, Rz, opt);
    h{j} = hj;
end

% compute fourier transform of every histogram
for j = 1 : Ns
    if use_gpu
        h{j} = fftn(gpuArray(h{j}));
    else
        h{j} = fftn(h{j});
    end
end

% compute cross-correlations
dx = zeros(Ns, Ns);
dy = zeros(Ns, Ns);
dz = zeros(Ns, Ns);
dm = false(Ns, Ns);
dx0 = zeros(Ns-1,1);
dy0 = zeros(Ns-1,1);
dz0 = zeros(Ns-1,1);

for i = 1 : Ns - 1
    hi = conj(h{i});
    for j = i + 1 : min(Ns, i + D) % either to Ns or maximally D more
        hj = ifftshift(real(ifftn(hi .* h{j})));
        xc = c(1);
        yc = c(2);
        zc = c(3);
        % centroid estimation
        gx = xc - 2*CR : xc + 2*CR;
        gy = yc - 2*CR : yc + 2*CR;
        gz = zc - 2*CR : zc + 2*CR;        
        d = gather(hj(gx, gy, gz));
        [gx, gy, gz] = ndgrid(gather(gx), gather(gy), gather(gz));
        gx = gx(:);
        gy = gy(:);
        gz = gz(:);
        d = d(:) - min(d(:));
        d = d / sum(d);
        for k = 1 : 20
            wc = exp(-4*log(2)*((xc - gx).^2+(yc-gy).^2+(zc-gz).^2)/CR^2);
            n = sum(wc .* d);
            xc = sum(gx .* d .* wc) / n;
            yc = sum(gy .* d .* wc) / n;
            zc = sum(gz .* d .* wc) / n;
        end
        sh = double([xc, yc, zc] - c);
        dx(i, j) = sh(1);
        dy(i, j) = sh(2);
        dz(i, j) = sh(3);
        dm(i, j) = true;
        if j == i + 1
            dx0(i) = sh(1);
            dy0(i) = sh(2);
            dz0(i) = sh(3);
        end
    end
end
[a,b] = ind2sub(size(dm), find(dm));
dx = dx(dm);
dy = dy(dm);
dz = dz(dm);
sx0 = cumsum(dx0);
sy0 = cumsum(dy0);
sz0 = cumsum(dz0);

% minimize cost function with some kind of regularization
options = optimoptions('fminunc', 'Display', 'none', 'MaxFunEvals', 1e5);
%options = optimset('Display', 'none', 'MaxFunEvals', 4e4, 'MaxIter', 4e4);

minimizer = @(x) lse_distance(x, a, b, dx, w, l);
sx = fminunc(minimizer, sx0, options);
sx = [0; sx];

minimizer = @(x) lse_distance(x, a, b, dy, w, l);
sy = fminunc(minimizer, sy0, options);
sy = [0; sy];

minimizer = @(x) lse_distance(x, a, b, dz, w, l);
sz = fminunc(minimizer, sz0, options);
sz = [0; sz];

% reduce by mean (so shift is minimal)
sx = sx - mean(sx);
sy = sy - mean(sy);
sz = sz - mean(sz);

% multiply by pixel size
sx = sx * sxyz;
sy = sy * sxyz;
sz = sz * sxyz;

% make gridded interpolants
Fx = griddedInterpolant(ti, sx, 'spline', 'linear');
Fy = griddedInterpolant(ti, sy, 'spline', 'linear');
Fz = griddedInterpolant(ti, sz, 'spline', 'linear');

% finally undrift
dx = Fx(t);
dy = Fy(t);
dz = Fz(t);

% also do it for every frame
ti = Rt(1):T/10:Rt(2);
dxt = Fx(ti);
dyt = Fy(ti);
dzt = Fz(ti);

fig = figure();
hold on;
plot(ti, dxt);
plot(ti, dyt);
plot(ti, dzt);
xlabel('time');
ylabel('drift');
legend('X', 'Y', 'Z');

end

function y = lse_distance(x, a, b, s, w, l)
% add shift 0 for the first frame
x = [0; x];
dx = x(b) - x(a);
y = w(b - a) .* (dx - s).^2;
y = mean(y);
% add regularization term (roughness)
r = l * sum(diff(x).^2);
y = y + r;
end