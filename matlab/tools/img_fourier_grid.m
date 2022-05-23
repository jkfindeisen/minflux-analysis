function [xi, yi, zi] = img_fourier_grid(dims, type_name)
% functions computed on this grid with center of mass at (0, 0) will if
% used in convolution via fft2 not produce any shift! so (0, 0) will be at
% the origin (1,1) and everything else as expected
%
% Syntax
%   [xi, yi] = img_fourier_grid(dims)
%
% 2012-07 by J.K.-F.

assert(nargin >= 1, 'Not enough arguments!');

if nargin < 2
    type_name = 'double';
end

number_dimensions = size(dims, 2);


switch number_dimensions
    case 1
        
        gx = cast(1 : dims(1), type_name);
        
        xi = ndgrid(gx);
        
        xi = ifftshift(xi);
        xi = xi - xi(1);
        
    case 2
        
        gx = cast(1 : dims(1), type_name);
        gy = cast(1 : dims(2), type_name);
        
        [xi, yi] = ndgrid(gx, gy);
        
        xi = ifftshift(xi);
        xi = xi - xi(1);
        
        yi = ifftshift(yi);
        yi = yi - yi(1);
        
    case 3
        
        gx = cast(1 : dims(1), type_name);
        gy = cast(1 : dims(2), type_name);
        gz = cast(1 : dims(3), type_name);
        
        [xi, yi, zi] = ndgrid(gx, gy, gz);
        
        xi = ifftshift(xi);
        xi = xi - xi(1);
        
        yi = ifftshift(yi);
        yi = yi - yi(1);
        
        zi = ifftshift(zi);
        zi = zi - zi(1);
        
    otherwise
        error('Unsupported dimensionality!');
end


end


% switch number_dimensions
%     case 2
%
%         % old code without ndgrid and ifftshift
%
%         nx = dims(1);
%         ny = dims(2);
%
%         if mod(nx, 2) == 0 % is even?
%             gridx = (-nx / 2) : 1 : (nx / 2 - 1);
%         else
%             gridx = -(nx - 1) / 2 : 1 : (nx - 1) / 2;
%         end
%
%         if mod(ny, 2) == 0 % is even?
%             gridy = (-ny / 2) : 1 : (ny / 2 - 1);
%         else
%             gridy = -(ny - 1) / 2 : 1 : (ny - 1) / 2;
%         end
%
%         [xi, yi] = ndgrid(gridx, gridy);
%         xi = circshift(xi, [gridx(1), 0]);
%         yi = circshift(yi, [0, gridy(1)]);
%         zi = [];
%
%     case 3
%
%         % here we use ndgrid and ifftshift
%         [xi, yi, zi] = ndgrid(1 : dims(1), 1 : dims(2), 1 : dims(3));
%         xi = ifftshift(xi);
%         xi = xi - xi(1, 1);
%         yi = ifftshift(yi);
%         yi = yi - yi(1, 1);
%         zi = ifftshift(zi);
%         zi = zi - zi(1, 1);
%
%     otherwise
%         error('Unsupported dimensionality!');
% end