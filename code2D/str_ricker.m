function [w,tw] = str_ricker(f,dt,tlength)
% str_RICKER: Ricker wavelet of central frequency f0.
%
% INPUT:
% f : central freq. in Hz (f <<1/(2dt) )
% dt: sampling interval in sec
% tlength : the duration of wavelet in sec
%
% OUTPUT: 
% w:  the Ricker wavelet
% tw: time axis
%
% Example
%
%   [w,tw] = str_ricker(10,0.004,0.2);
%    plot(tw,w);
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

if nargin==3
    nw=floor(tlength/dt)+1;
else
    nw=2.2/f/dt;
    nw=2*floor(nw/2)+1;
end
nc=floor(nw/2);
w = zeros(nw,1);

k=[1:1:nw]';

alpha = (nc-k+1).*f*dt*pi;
beta=alpha.^2;
w = (1.-beta.*2).*exp(-beta);

if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
end