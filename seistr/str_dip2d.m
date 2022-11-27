function [dip] = str_dip2d(din,niter,liter,order,eps_dv, eps_cg, tol_cg,rect,verb)
% str_dip2d: dip estimation based on shaping regularized PWD algorithm 
% (independent implementation)
%
% INPUT
% din: input data (nt*nx)
% niter: number of nonlinear iterations
% liter: number of linear iterations (in divn)
% order: accuracy order
% eps_dv: eps for divn  (default: 0.01)
% eps_cg: eps for CG    (default: 1)
% tol_cg: tolerence for CG (default: 0.000001)
% rect:  smoothing radius (ndim*1)
% verb: verbosity flag
%
% OUTPUT
% dip:  2D slope
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

if nargin==1
    niter = 5;
    liter = 20;
    eps_dv = 0.01;
    eps_cg = 1;
    tol_cg = 0.000001;
    order=2;
    rect(1) = 10;
    rect(2) = 10;
    rect(3) = 1;
    verb=1;
end

p0 = 0.0;
nj1 = 1;
nj2 = 1;

dim = 3;

n = zeros(dim, 1);
n1 = size(din,1);
n2 = size(din,2);
n(1) = n1;
n(2) = n2;
n(3) = 1;

n123 = n(1) * n(2) * n(3);

ratio=zeros(size(din));
dip=zeros(size(din));

for iter=1:niter
%corresponding to the eq.21 in the paper 
[u1,u2] = str_conv_allpass(din,dip,order); 
[ ratio ] = str_divne(-u2, u1, liter, rect, n, eps_dv, eps_cg, tol_cg,verb);
dip=dip+ratio;
end

return



