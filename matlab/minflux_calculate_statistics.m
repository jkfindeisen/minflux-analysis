function stat = minflux_calculate_statistics(minflux)
% Calculates additional statistics on a gathered Minflux experiment
%
% All physical lengths are in m.

assert(nargin == 1);

%% extract all relevant information from the minflux data structure, remove all invalid entries
valid = minflux.vld(:);
stat.fraction_valid = mean(valid);
efc = minflux.itr.efc;
efo = minflux.itr.efo;
fbg = minflux.itr.fbg(:, end);
pos = squeeze(minflux.itr.loc(:, end, :));
t = minflux.tim(:);
t = t - t(1); % w.r.t. start at zero
id = minflux.tid(:);

%% center frequency ratio

% average NaN, Inf rate of columns in EFC, EFO
r1 = mean(isnan(efc) | isinf(efc));
r2 = mean(isnan(efo) | isinf(efo));
r3 = r1 < 0.5 & r2 < 0.5;
idx = find(r3, 1, 'last');
if idx < size(efc, 2)
    fprintf('take cfr from iteration %d\n', idx);
end
cfr = efc(:, idx) ./ efo(:, idx); % CFR is now holding the idx iteration

%% ranges
x = pos(:, 1); % in m
y = pos(:, 2);
z = pos(:, 3);

% reasonable boundaries (will also be used further below)
a = 0.01;
Rx = quantile(x, [a, 1-a]);
Ry = quantile(y, [a, 1-a]);
Rz = quantile(z, [a, 1-a]);
% minimal boundaries 200nm (because of CR in drift_correction)
Rmin = 200e-9;
if diff(Rx) < Rmin
    Rx = Rx + (Rmin - diff(Rx))/2*[-1,1];
end
if diff(Ry) < Rmin
    Ry = Ry + (Rmin - diff(Ry))/2*[-1,1];
end
if diff(Rz) < Rmin && diff(Rz) > 1e-9 % only if 3D
    Rz = Rz + (Rmin - diff(Rz))/2*[-1,1];
end
Rt = t([1,end]);

%% drift correction
T = numel(unique(id))*diff(Rx)*diff(Ry)/3e6; % heuristic for optimal length of time window
T = min([T, diff(Rt)/2, 3600]); % need at least two time windows
T = max([T, 600]); % but at least 10 minutes long
sxy = 2e-9;
sxyz = 5e-9;
use_gpu = true;

% switch for 2D/3D
is3D = diff(Rz) > 1e-9;

if is3D
    % 3D drift correction
    [dx, dy, dz, dxt, dyt, dzt, ti, fig] = drift_correction_time_windows_3D(x, y, z, t, Rx, Ry, Rz, T, sxyz, use_gpu);
    close(fig);
    dpos = pos - [dx, dy, dz];
    drift = [ti(:), dxt(:), dyt(:), dzt(:)];
else
    % 2D drift correction
    [dx, dy, dxt, dyt, ti, fig] = drift_correction_time_windows_2D(x, y, t, Rx, Ry, T, sxy, use_gpu);
    close(fig);
    dpos = pos - [dx, dy, zeros(size(dx))];
    drift = [ti(:), dxt(:), dyt(:)];
end

% will work on dpos from now one

%% combine localization of each event, compute std deviation in x,y,z
[~, ~, uid] = unique(id);
c.n = accumarray(uid, 1);
N = numel(c.n);
c.t = accumarray(uid, t) ./ c.n;
c.pos = zeros(N, 3);
c.std_xyz = zeros(N, 3);
for i = 1 : 3
    c.pos(:, i) = accumarray(uid, dpos(:, i)) ./ c.n;
    c.std_xyz(:, i) = accumarray(uid, dpos(:, i), [], @std);
end

% compute std-error in x,y and in z
Tn = 3;
s = c.std_xyz(c.n >= Tn, :); % ignore only seen once (likely false positives anyway)
n = c.n(c.n >= Tn);
c.ste_xyz = s ./ sqrt(n); % standard error = standard deviation / sqrt(N)

%% start and end times of each binding event (first and last time)
stat.t_start = accumarray(uid, t, [], @min);
stat.t_end = accumarray(uid, t, [], @max);
stat.t_scan = stat.t_start(2:end, 1) - stat.t_end(1:end-1); % start of next - end of previous
stat.t_loc = stat.t_end - stat.t_start; % end - start

%% fourier ring correlation (only in x,y)
R = 8; % five repetitions to get more stable results

% x, y coordinates (in nm)
x = dpos(:, 1);
y = dpos(:, 2);
for r = 1 : R
    % new split
    ix = rand(size(x)) < 0.5;
    
    % 2D render of x,y (histogram, 1nm pixel size)
    sxy = 1e-9;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    [estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy);
    frc.resolution(r) = estimated_resolution;
    frc.qi{r} = qi;
    frc.ci{r} = ci;
    
    % now on combined
    x = c.pos(:, 1);
    y = c.pos(:, 2);
    ix = rand(size(x)) < 0.5;
    h1 = render_xy(x(ix), y(ix), sxy, sxy, Rx, Ry);
    h2 = render_xy(x(~ix), y(~ix), sxy, sxy, Rx, Ry);
    
    [estimated_resolution, ~, qi, ci] = img_fourier_ring_correlation(h1, h2, size(h1)*sxy);
    frc_combined.resolution(r) = estimated_resolution;
    frc_combined.qi{r} = qi;
    frc_combined.ci{r} = ci;
end

%% assemble output
stat.pos = pos;
stat.dpos = dpos; % drift corrected position
stat.drift = drift;
stat.combined = c;
stat.cfr = cfr;
stat.t = t;
stat.id = id;
stat.fbg = fbg;
stat.is3D = is3D;
stat.frc = frc;
stat.frc_combined = frc_combined;

end