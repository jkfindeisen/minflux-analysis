function [resolution, fc] = img_fourier_correlation_3d(image1, image2, physical_lengths, kernel_size)
% Calculates the Fourier correlation and then computes an ellipsoid at the
% 3D Fourier correlation. Expects the physical image length in meter.

assert(nargin == 4);

% compute smoothing kernel
L = ceil(3 * kernel_size);
[x, y, z] = ndgrid(-L(1) : L(1),-L(2) : L(2), -L(3) : L(3));
k = exp(-log(2) * 4 * (x.^2 / kernel_size(1)^2 + y.^2 / kernel_size(2)^2 + z.^2 / kernel_size(3)^2));
fc_smoothing_kernel = k / sum(k(:));
clear x y z k L kernel_size;

% fourier transforms
f1 = fftshift(fftn(image1));
f2 = fftshift(fftn(image2));
clear image1 image2;

% derived quantities
a = f1 .* conj(f2);
b = f1 .* conj(f1);
c = f2 .* conj(f2);
clear f1 f2;

% convolve and then compute correlation
a_sm = convn(a, fc_smoothing_kernel, 'same');
b_sm = convn(b, fc_smoothing_kernel, 'same');
c_sm = convn(c, fc_smoothing_kernel, 'same');
clear a b c;
fc = a_sm ./ sqrt(b_sm .* c_sm);
fc = real(fc);
clear a_sm b_sm c_sm;

% should be separable, try addup in other two directions each time
fcx = sum(sum(fc, 2), 3);
fcx = fcx(:);
nx = length(fcx);

fcy = sum(sum(fc, 1), 3);
fcy = fcy(:);
ny = length(fcy);

fcz = sum(sum(fc, 1), 2);
fcz = fcz(:);
nz = length(fcz);

% get center of mass
cx = round(sum(fcx .* (1 : nx).') / sum(fcx));
cy = round(sum(fcy .* (1 : ny).') / sum(fcy));
cz = round(sum(fcz .* (1 : nz).') / sum(fcz));

% fold together
lx = min(cx, nx - cx + 1)-1;
ly = min(cy, ny - cy + 1)-1;
lz = min(cz, nz - cz + 1)-1;

fx = flip(fcx(cx-lx:cx)) + fcx(cx:cx+lx);
fy = flip(fcy(cy-ly:cy)) + fcy(cy:cy+ly);
fz = flip(fcz(cz-lz:cz)) + fcz(cz:cz+lz);

% smooth fx, fy
fxs = sgolayfilt(fx, 1, 27); % TODO adapt this on the real length scale
fys = sgolayfilt(fy, 1, 27);

% normalize
fxn = fxs / fxs(1);
fyn = fys / fys(1);
fzn = fz / fz(1);

% find threshold 1/7
t = 1 / 7;
tx = omex_calcHLT(fxn, false, t);
ty = omex_calcHLT(fyn, false, t);
tz = omex_calcHLT(fzn, false, t);

% convert this to resolution
resolution = physical_lengths ./ [tx, ty, tz];

end