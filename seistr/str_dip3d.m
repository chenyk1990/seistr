function [dip_i,dip_x] = str_dip3d(din,niter,liter,order,eps_dv, eps_cg, tol_cg,rect,verb)
% str_dip3d: 3D dip estimation based on shaping regularized PWD algorithm
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
% dipi:  inline 3D slope
% dipx:  xline 3D slope
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

if nargin==1
    niter = 5;
    liter = 10;
    eps_dv = 0.01;
    eps_cg = 1;
    tol_cg = 0.000001;
    order=2;
    rect(1) = 5;
    rect(2) = 5;
    rect(3) = 5;
    verb=1;
end

if nargin==3
    eps_dv = 0.01;
    eps_cg = 1;
    tol_cg = 0.000001;
    order=2;
    rect(1) = 5;
    rect(2) = 5;
    rect(3) = 5;
    verb=1;
end

dim = 3;
n = zeros(dim, 1);
n1 = size(din,1);
n2 = size(din,2);
n3 = size(din,3);
n(1) = n1;
n(2) = n2;
n(3) = n3;

n123 = n(1) * n(2) * n(3);


dip_i=zeros(size(din));
dip_x=zeros(size(din));

for iter=1:niter
    
    %corresponding to the eq.21 in the paper
    [u1_i,u2_i] = str_conv_allpass_i(din,dip_i,order); % inline linearization using the updated dip    
    [ ratio_i ] = str_divne(-u2_i, u1_i, liter, rect, n, eps_dv, eps_cg, tol_cg,verb);

    %corresponding to the eq.21 in the paper
    [u1_x,u2_x] = str_conv_allpass_x(din,dip_x,order); % xline linearization using the updated dip
    [ ratio_x ] = str_divne(-u2_x, u1_x, liter, rect, n, eps_dv, eps_cg, tol_cg,verb);
    
    dip_i=dip_i+ratio_i;
    dip_x=dip_x+ratio_x;
end

return