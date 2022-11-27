function ds = str_pwsmooth_lop2d(dn,dip,ns,order,eps)
% str_pwsmooth_lop: plane-wave smoothing 
%
% INPUT:
% dn: model   noisy data
% dip: slope (2D array)
% ns:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);

% OUTPUT:
% ds:  smoothed data
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

n1=size(dn,1);
n2=size(dn,2);

ns2=2*ns+1;%spray diameter

utmp=str_pwspray_lop2d(dn,dip,ns,order,eps);

u=reshape(utmp,n1,ns2,n2);

ds=squeeze(sum(u,2)/ns2);

return




