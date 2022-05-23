function [estimated_resolution, fc, qi, ci] = img_fourier_ring_correlation(image1, image2, physical_image_size, frc_smoothing_kernel)
% Does fourier ring correlation analysis on STED images via splitting the
% data in one direction.
%
% Expects physical_image_size in [m]!!

assert(nargin >= 3, 'Not enough arguments!');

if nargin < 4 || isempty(frc_smoothing_kernel)
    frc_smoothing_kernel = fspecial('gaussian', [31, 31], 1);
end

% fourier transform
f1 = fft2(image1);
f2 = fft2(image2);

% derived quantities for correlation
a = f1 .* conj(f2);
b = f1 .* conj(f1);
c = f2 .* conj(f2);

% 2D image representation (first smooth, then real/absolute value)
a_sm = conv2(ifftshift(a), frc_smoothing_kernel, 'same');
b_sm = conv2(ifftshift(b), frc_smoothing_kernel, 'same');
c_sm = conv2(ifftshift(c), frc_smoothing_kernel, 'same');
fc = a_sm ./ sqrt(b_sm .* c_sm);
fc = real(fc);

% calculate frequency space grid
[qx, qy] = img_fourier_grid(size(image1));
qx = qx / physical_image_size(1);
qy = qy / physical_image_size(2);
q = sqrt(qx.^2 + qy.^2);

% bin a,b,c in dependence of q
B = 5e5; % bin size (in pixel in fourier space) % TODO must be adapted to physical stack size (this is for m)!!
qi = round(q / B);
idx = qi(:) + 1;
qi = (0 : max(qi(:))).' * B;
aj = accumarray(idx, a(:));
bj = accumarray(idx, b(:));
cj = accumarray(idx, c(:));
ci = double(real(aj ./ sqrt(bj .* cj)));
idx = qi < max(qi)*0.8; % cut a bit
qi = qi(idx);
ci = ci(idx);

% this is additional smoothing
if license('test', 'Signal_Toolbox') % not everyone has it
    ci = sgolayfilt(ci, 1, 7);
else
    ci = conv(ci, [1, 1, 1] / 3, 'same');
end

% determine resolution
q_critical = qi(find(ci < 1 / 7, 1, 'first'));
if isempty(q_critical)
    q_critical = qi(end);
end
estimated_resolution = 1 / q_critical;

end