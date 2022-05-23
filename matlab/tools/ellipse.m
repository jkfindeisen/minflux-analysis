function f = ellipse(x, y, p)
% elliptic ring with 6 parameters
%   p(1) - center x
%   p(2) - center y
%   p(3) - angle of first semi axis to x axis [rad]
%   p(4) - first semi axis length
%   p(5) - second semi axis length
%   p(6) - width of elliptic ring

assert(nargin == 3);

% move coordinate system (COS) to center point
xc = x - p(1);
yc = y - p(2);

% rotate COS with rotation angle
cp = cos(p(3));
sp = sin(p(3));
xr = cp * xc - sp * yc;
yr = sp * xc + cp * yc;

% scale directions by semi axes
xr = xr / p(4);
yr = yr / p(5);

% normalized radius
r = sqrt(xr.^2 + yr.^2);

% distance towards ellipse
d = (r - 1).^2;

% gaussian function
f = exp(-4*log(2)*d/(p(6)/sqrt(p(4)*p(5)))^2);

% normalization
f = f / sum(f(:));

end