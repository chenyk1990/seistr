function [u1,u2] = str_conv_allpass(din,dip,order)
% str_conv_allpass: Convolutional operator implemented by an allpass filter
% 
% Linearized inverse problem
% C'(\sigma)d\Delta \sigma =  C(\sigma)d
% 
% IPNUT:
% din: input data
% dip: 3D dip
% order: accuracy order
%
% OUTPUT:
% u1: C'(\sigma)d (denominator)
% u2: C(\sigma)d  (numerator)
% 
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

u1=zeros(size(din));
u2=zeros(size(din));

[n1,n2]=size(din);

if order==1
    nw=1;
else
    nw=2;
end
filt1=zeros(2*nw+1,1);
filt2=zeros(2*nw+1,1);

% cprresponding to eqs.17-19
% destruct the data
for i1=nw+1:n1-nw
    for i2=1:n2-1
        
        if order==1
            filt1=B3d(dip(i1,i2));
            filt2=B3(dip(i1,i2));
        else
            filt1=B5d(dip(i1,i2));
            filt2=B5(dip(i1,i2));
        end
        for iw=-nw:nw
            u1(i1,i2)=u1(i1,i2)+(din(i1+iw,i2+1)-din(i1-iw,i2))*filt1(iw+nw+1);
            u2(i1,i2)=u2(i1,i2)+(din(i1+iw,i2+1)-din(i1-iw,i2))*filt2(iw+nw+1);
        end
    end
end

return

% cprresponding to eqs.13-16
% form the filters
function [b3 ] = B3(sigma)
% B3 coefficient
% sigma: slope

b3(1)=(1-sigma)*(2-sigma)/12;
b3(2)=(2+sigma)*(2-sigma)/6;
b3(3)=(1+sigma)*(2+sigma)/12;

return

function [b3d ] = B3d(sigma)
% B3 coefficient derivative
% sigma: slope
b3d(1)=-(2-sigma)/12-(1-sigma)/12;
b3d(2)=(2-sigma)/6-(2+sigma)/6;
b3d(3)=(2+sigma)/12+(1+sigma)/12;

return

function [b5 ] = B5(sigma)
% B5 coefficient
% sigma: slope

b5(1)=(1-sigma)*(2-sigma)*(3-sigma)*(4-sigma)/1680;
b5(2)=(4-sigma)*(2-sigma)*(3-sigma)*(4+sigma)/420;
b5(3)=(4-sigma)*(3-sigma)*(3+sigma)*(4+sigma)/280;
b5(4)=(4-sigma)*(2+sigma)*(3+sigma)*(4+sigma)/420;
b5(5)=(1+sigma)*(2+sigma)*(3+sigma)*(4+sigma)/1680;

return

function [b5d ] = B5d(sigma)
% B5 coefficient derivative
% sigma: slope

b5d(1)=-(2-sigma)*(3-sigma)*(4-sigma)/1680-...
    (1-sigma)*(3-sigma)*(4-sigma)/1680-...
    (1-sigma)*(2-sigma)*(4-sigma)/1680-...
    (1-sigma)*(2-sigma)*(3-sigma)/1680;

b5d(2)=-(2-sigma)*(3-sigma)*(4+sigma)/420-...
    (4-sigma)*(3-sigma)*(4+sigma)/420-...
    (4-sigma)*(2-sigma)*(4+sigma)/420+...
    (4-sigma)*(2-sigma)*(3-sigma)/420;

b5d(3)=-(3-sigma)*(3+sigma)*(4+sigma)/280-...
    (4-sigma)*(3+sigma)*(4+sigma)/280+...
    (4-sigma)*(3-sigma)*(4+sigma)/280+...
    (4-sigma)*(3-sigma)*(3+sigma)/280;

b5d(4)=-(2+sigma)*(3+sigma)*(4+sigma)/420+...
    (4-sigma)*(3+sigma)*(4+sigma)/420+...
    (4-sigma)*(2+sigma)*(4+sigma)/420+...
    (4-sigma)*(2+sigma)*(3+sigma)/420;

b5d(5)=(2+sigma)*(3+sigma)*(4+sigma)/1680+...
    (1+sigma)*(3+sigma)*(4+sigma)/1680+...
    (1+sigma)*(2+sigma)*(4+sigma)/1680+...
    (1+sigma)*(2+sigma)*(3+sigma)/1680;

return