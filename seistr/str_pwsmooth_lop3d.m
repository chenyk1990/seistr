function ds=str_pwsmooth_lop3d(dn,dipi,dipx,r1,r2,eps,order)
% str_pwsmooth_lop3d: 3D plane-wave smoothing
%
% INPUT:
% dn: model  noisy data
% dipi: inline slope
% dipx: xline slope
% r1,r2:    spray radius
% order:    PWD order
% eps: regularization (default:0.01);
% 
% OUTPUT:
% ds: smoothed data
% 
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

np=(2*r1+1)*(2*r2+1);
%flattening
[u] = str_pwspray_lop3d(dn,dipi,dipx,r1,r2,order,eps);

%smoothing
ds=squeeze(sum(u,2)/np);






