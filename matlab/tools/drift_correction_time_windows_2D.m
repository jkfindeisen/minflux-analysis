function [dx, dy, dxt, dyt, ti, fig] = drift_correction_time_windows_2D(x, y, t, Rx, Ry, T, sxy, use_gpu)
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
CR = 10; % center of mass
D = Ns; % how many frames maximally in the cross-correlation
w = linspace(1, 0.5, D).'; % weight with distance
l = 1e-1;

% get dimensions
Nx = ceil(diff(Rx) / sxy);
Ny = ceil(diff(Ry) / sxy);
c = round([Nx, Ny]/2) + 1;

% create all the histograms
h = cell(Ns, 1);
ti = zeros(Ns, 1); % average times of the snapshots
opt = struct('type', 'fixed_gaussian', 'fwhm', 3*sxy);
for j = 1 : Ns
    t0 = Rt(1) + (j-1)*T;
    idx = t >= t0 & t < t0 + T;
    ti(j) = mean(t(idx));
    hj = render_xy(x(idx), y(idx), sxy, sxy, Rx, Ry, opt);
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
dm = false(Ns, Ns);
dx0 = zeros(Ns-1,1);
dy0 = zeros(Ns-1,1);

for i = 1 : Ns - 1
    hi = conj(h{i});
    for j = i + 1 : min(Ns, i + D) % either to Ns or maximally D more
        hj = ifftshift(real(ifftn(hi .* h{j})));
        xc = c(1);
        yc = c(2);
        % centroid estimation
        gx = xc - 2*CR : xc + 2*CR;
        gy = yc - 2*CR : yc + 2*CR;
        d = gather(hj(gx, gy));
        [gx, gy] = ndgrid(gather(gx), gather(gy));
        gx = gx(:);
        gy = gy(:);
        d = d(:) - min(d(:));
        d = d / sum(d);
        for k = 1 : 20
            wc = exp(-4*log(2)*((xc - gx).^2+(yc-gy).^2)/CR^2);
            n = sum(wc .* d);
            xc = sum(gx .* d .* wc) / n;
            yc = sum(gy .* d .* wc) / n;
        end
        sh = double([xc, yc] - c);
        dx(i, j) = sh(1);
        dy(i, j) = sh(2);
        dm(i, j) = true;
        if j == i + 1
            dx0(i) = sh(1);
            dy0(i) = sh(2);
        end
    end
end
[a,b] = ind2sub(size(dm), find(dm));
dx = dx(dm);
dy = dy(dm);
sx0 = cumsum(dx0);
sy0 = cumsum(dy0);

% minimize cost function with some kind of regularization
options = optimoptions('fminunc', 'Display', 'none', 'MaxFunEvals', 1e5);
%options = optimset('Display', 'none', 'MaxFunEvals', 4e4, 'MaxIter', 4e4);

minimizer = @(x) lse_distance(x, a, b, dx, w, l);
sx = fminunc(minimizer, sx0, options);
sx = [0; sx];

minimizer = @(x) lse_distance(x, a, b, dy, w, l);
sy = fminunc(minimizer, sy0, options);
sy = [0; sy];

% reduce by mean (so shift is minimal)
sx = sx - mean(sx);
sy = sy - mean(sy);

% multiply by pixel size
sx = sx * sxy;
sy = sy * sxy;

% make gridded interpolants
Fx = griddedInterpolant(ti, sx, 'spline', 'linear');
Fy = griddedInterpolant(ti, sy, 'spline', 'linear');

% finally undrift
dx = Fx(t);
dy = Fy(t);

% also do it for every frame
ti = Rt(1):T/10:Rt(2);
dxt = Fx(ti);
dyt = Fy(ti);

fig = figure();
hold on;
plot(ti, dxt);
plot(ti, dyt);
xlabel('time');
ylabel('drift');
legend('X', 'Y');

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