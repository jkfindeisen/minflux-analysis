function [fc] = img_fourier_correlation(image1, image2, kernel_size)
% Does fourier ring correlation analysis on STED images via splitting the
% data in one direction.
%
% Expects physical_image_size in [m]!!

assert(nargin == 3, 'Not enough arguments!');

L = ceil(3 * kernel_size);
if length(kernel_size) == 1
    [x, y] = ndgrid(-L:L,-L:L);
    k = exp(-(x.^2+y.^2) / (2 * kernel_size^2));
else
    [x, y] = ndgrid(-L(1):L(1),-L(2):L(2));
    k = exp(-(x.^2 / kernel_size(1)^2 + y.^2 / kernel_size(2)^2) / 2);
end
fc_smoothing_kernel = k / sum(k(:));
% fc_smoothing_kernel = fspecial('gaussian', ceil(8 * kernel_size) * [1,1] , kernel_size);

%image1 = image1 - mean(image1(:));
%image2 = image2 - mean(image2(:));

% fourier transform
f1 = fft2(image1);
f2 = fft2(image2);

% derived quantities for correlation
a = f1 .* conj(f2);
b = f1 .* conj(f1);
c = f2 .* conj(f2);

% 2D image representation (first smooth, then real/absolute value)
a_sm = conv2(ifftshift(a), fc_smoothing_kernel, 'same');
b_sm = conv2(ifftshift(b), fc_smoothing_kernel, 'same');
c_sm = conv2(ifftshift(c), fc_smoothing_kernel, 'same');
fc = a_sm ./ sqrt(b_sm .* c_sm);
fc = real(fc);

end