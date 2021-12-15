function [dout] = str_pwspray_lop3d(din,dipi,dipx,ns2,ns3,order,eps)
%str_pwspray_lop3d: 3D plane-wave spray operator 
% 
% INPUT:
% din: input
% dipi: inline slope
% dipx: xline slope
% ns2/3: smoothing radius
% order: PWD order
% eps: regularization (default:0.01);
% 
% OUTPUT:
% dout: result
% 
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.



[n1,n2,n3]=size(din);
n23=n2*n3;


np2=2*ns2+1;%spray diameter
np3=2*ns3+1;%spray diameter
np=np2*np3;

trace=zeros(n1,1);
e=eps*eps;
nw=order;

% traveltime table
t=zeros(np2,np3);
for i3=0:np3-1
    for i2=0:np2-1
	    t(i2+1,i3+1) = hypot(i2-ns2,i3-ns3); %square root of the sum of squares
    end
end

[~,visit]=sort(t(:));

u=zeros(n1,np,n23);
din=reshape(din,n1,n23);
p1=reshape(dipi,n1,n23);
p2=reshape(dipx,n1,n23);

for ii=0:n23-1
    
    
    u(:,ns3*np2+ns2+1,ii+1)=din(:,ii+1); %central trace in the predicted dimension
    
    i2=mod(ii,n2);
    i3=floor(ii/n2);
    
    for ip=0:np-1

        [up, jp, up2, up3 ] = get_update(ip, np2, np3, t, visit);
        
        % from jp to j
        k2=mod(jp,np2);
        k3=floor(jp/np2);
        
        j2=i2+k2-ns2;
        j3=i3+k3-ns3;
                
        if j2<0 || j2>=n2 || j3<0 || j3>=n3
            continue;
        end
        
        j=j2+j3*n2;
        
        if bitand(up,1)
            if up2
                if j2==0
                    continue;
                end
                j2=j-1;
                q2=p1(:,j2+1);
                k2=jp-1;
            else
                if j2==n2-1
                    continue;
                end
                j2=j+1;
                q2=p1(:,j+1);
                k2=jp+1;
            end
        end
        
        if bitand(up,2)
            if up3
                if j3==0
                    continue;
                end
                j3=j-n2;
                q3=p2(:,j3+1);
                k3=jp-np2;
            else
                if j3==n3-1
                    continue;
                end
                j3=j+n2;
                q3=p2(:,j+1);
                k3=jp+np2;
            end
        end
        

        %corresponding to the eqs.24/25
        %prediction forward or backward
        switch up
            
            case 0
                continue;
            case 1

                [u(:,jp+1,j+1)] = predict1_step(e,nw,up2,q2,u(:,k2+1,j2+1));
            case 2
                [u(:,jp+1,j+1)] = predict1_step(e,nw,up3,q3,u(:,k3+1,j3+1));
            case 3
                [u(:,jp+1,j+1)] = predict2_step(e,nw,up2,up3,q2,q3,u(:,k2+1,j2+1),u(:,k3+1,j3+1));
        end
        
        
    end
end
dout=reshape(u,n1,np,n2,n3);

return

function [up, jj, up1, up2 ] = get_update(ii, n1,n2, t, visit)
% forward or backward options
% 
% INPUT
% ii: loop number
% n1:   [n1,n2]=size(t)
% n2:   [n1,n2]=size(t)
% t:    traveltime
% visit: sorted sequence of t
% 
% OUTPUT
% jj: local linear index
% up: prediction type
% up1: forward or backward
% up2: forward or backward
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.


n12=n1*n2;


jj=visit(ii+1)-1;
t1=t(jj+1);

i1=mod(jj,n1);
i2=floor(jj/n1);

up=0;
if n1>1
    a1=jj-1;
    b1=jj+1;
    
    up1= (i1 && (i1==n1-1 || 1~= (t(a1+1)>t(b1+1))));
    if up1
        c1=a1;
    else
        c1=b1;
    end
    if t1>t(c1+1)
        up = bitor(up, 1); %Bit-wise OR.
    end
end

if n2>1
    a2=jj-n1;
    b2=jj+n1;
    
    up2= (i2 && (i2==n2-1 || 1~= (t(a2+1)>t(b2+1))));
    if up2
        c2=a2;
    else
        c2=b2;
    end
    if t1>t(c2+1)
        up = bitor(up, 2); %Bit-wise OR.
    end
end

return

function [trace] = predict1_step(e,nw,forw,pp,trace1)
% predict1_step: prediction step from one trace
% 
% INPUT:
% e: regularization parameter (default, 0.01*0.01);
% nw: accuracy order
% forw: forward or backward
% pp: slope
% trace1: input trace
% 
% OUTPUT:
% trace: output trace
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.


n1=size(trace1,1);
nb=2*nw;
eps=e;
eps2=e;

diag=zeros(n1,1);
offd=zeros(n1,nb);


[diag,offd] = regularization(diag,offd,nw,eps,eps2);



[w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp);


t0=trace1(1);
t1=trace1(2);
t2=trace1(n1-1);
t3=trace1(n1);

[trace] = pwd_set(0,w,diag,offd,pp,trace1);

trace(1)=trace(1)+eps2*t0;
trace(2)=trace(2)+eps2*t1;

trace(n1-1)=trace(n1-1)+eps2*t2;
trace(n1)=trace(n1)+eps2*t3;


trace = banded_solve(n1,nb,diag,offd,trace);   


return


function [diag,offd] = regularization(diag,offd,nw,eps,eps2)
% fill diag and offd using regularization 
% 
% INPUT:
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% eps: regularization parameter (default: e*e, e=0.01);
% eps2: second regularization parameter (default, same as eps)
% nw: accuracy order (nb=2*nw)
%
% OUTPUT:
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

%
nb=2*nw;
n1=length(diag);

for i1=0:n1-1
    diag(i1+1)=6*eps;
    offd(i1+1,1)=-4*eps;
    offd(i1+1,2)=eps;
    
    for ib=2:nb-1
        offd(i1+1,ib+1)=0.0;
    end
end

    diag(1)=eps2+eps;
    diag(2)=eps2+5*eps;
    
    diag(n1)=eps2+eps;
    diag(n1-1)=eps2+5*eps;  
    
    offd(1,1)=-2*eps;
    offd(n1-1,1)=-2*eps;

return

function [w,diag,offd] = pwd_define(forw,diag,offd,n1,nw,pp)
% pwd_define: matrix multiplication
% 
% INPUT:
% forw:forward or backward
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% n1: trace length
% nw: PWD filter(accuracy) order (default nw=1)
% pp: slope                  (1D array)

% OUTPUT:
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% w: PWD object(struct)
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.


%define PWD object(struct)
w=struct;
w.n=n1;
w.na=2*nw+1;
w.a=zeros(n1,w.na);%2D array
w.b=zeros(w.na,1);%1D array


n=w.n;
nb=2*nw;


for i=0:n-1
    w.b=passfilter(pp(i+1),nw);
    for j=0:w.na-1
        if(forw)
            w.a(i+1,j+1)=w.b(w.na-j);
        else
            w.a(i+1,j+1)=w.b(j+1);
        end
        
    end
end

for i=0:n-1
    for j=0:w.na-1
        k=i+j-nw;
        if k>=nw && k<n-nw
            aj=w.a(k+1,j+1);
            diag(i+1)=diag(i+1)+aj*aj;
        end
    end
    for m=0:2*nw-1
        for j=m+1:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                aj=w.a(k+1,j+1);
                am=w.a(k+1,j-m);
                offd(i+1,m+1)=offd(i+1,m+1)+am*aj;
            end
        end
    end
    
end

return


function [a,b] = passfilter(p,nw)
% passfilter: find filter coefficients 
% All-pass plane-wave destruction filter coefficients
% 
% INPUT:
% p: slope
% nw: PWD filter(accuracy) order
%
% OUTPUT:
% a: output filter (n+1) (1D array)
% b: temp variable of a
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.


n=nw*2;
b=zeros(n+1,1);
for k=0:n
    bk=1;
    for j=0:n-1
        if (j<n-k)
            bk=bk*(k+j+1.0)/(2*(2*j+1)*(j+1));
        else
            bk=bk*1.0/(2*(2*j+1));
        end
        
    end
    b(k+1)=bk;
end

for k=0:n
    ak=b(k+1);
    for j=0:n-1
        if j<n-k
            ak=ak*(n-j-p);
        else
            ak=ak*(p+j+1);
        end
    end
    a(k+1)=ak;
end
return


function [out] = pwd_set(adj,w,diag,offd,pp,inp)
% pwd_set: matrix multiplication
% 
% INPUT:
% adj:adjoint flag
% w: PWD object(struct)
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% pp: slope                  (1D array)
% inp: model
%
% OUTPUT:
% out: data 
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.




n=w.n;
nw=(w.na-1)/2;

%% pwd_set
tmp=zeros(n,1);

if adj
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=0:n-1
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                tmp(k+1)=tmp(k+1)+w.a(k+1,j+1)*inp(i+1);
            end
            
        end
    end
    
    for i=0:n-1
        out(i+1)=0.0;
    end
    
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            out(k+1)=out(k+1)+w.a(i+1,j+1)*tmp(i+1);
        end
    end
else
    for i=0:n-1
        tmp(i+1)=0.0;
    end
    for i=nw:n-nw-1
        for j=0:w.na-1
            k=i+j-nw;
            tmp(i+1)=tmp(i+1)+w.a(i+1,j+1)*inp(k+1);
        end
    end
    for i=0:n-1
        out(i+1)=0.0;
        for j=0:w.na-1
            k=i+j-nw;
            if k>=nw && k<n-nw
                out(i+1)=out(i+1)+w.a(k+1,j+1)*tmp(k+1);
            end
        end 
    end
    
end


return


function b = banded_solve(n,band,diag,offd,b)
% banded_solve: Banded matrix solver
% 
% INPUT:
% n:    matrix size
% band: band size
% diag: defined diagonal .   (1D array)
% offd: defined off-diagonal (2D array)
% b: input trace
%
% OUTPUT:
% b: trace solution
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.
%


%define Band object(struct)
slv=struct;
slv.n=n;
slv.band=band;
slv.d=zeros(n,1);       %1D array
slv.o=zeros(n-1,band);  %2D array


% define the banded matrix
for k=0:slv.n-1
    t=diag(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*slv.o(k-m,m+1)*slv.d(k-m);
    end
    slv.d(k+1)=t;
    n1=min(slv.n-k-1,slv.band);
    for n=0:n1-1
        t=offd(k+1,n+1);
        m1=min(k,slv.band-n-1);
        for m=0:m1-1
            t=t-slv.o(k-m,m+1)*slv.o(k-m,n+m+2)*slv.d(k-m);
        end
        slv.o(k+1,n+1)=t/slv.d(k+1);
    end
    
    
end

% the solver 
for k=1:slv.n-1
    
    t=b(k+1);
    m1=min(k,slv.band);
    for m=0:m1-1
        t=t-slv.o(k-m,m+1)*b(k-m);
    end
    b(k+1)=t;
end
for k=slv.n-1:-1:0
    t=b(k+1)/slv.d(k+1);
    m1=min(slv.n-k-1,slv.band);
    for m=0:m1-1
        t=t-slv.o(k+1,m+1)*b(k+m+2);
    end
    b(k+1)=t;
end



return


function [trace] = predict2_step(e,nw,forw1,forw2,pp1,pp2,trace1,trace2)
% predict2_step: prediction step from two trace
% 
% INPUT:
% e: regularization parameter (default, 0.01*0.01);
% nw: accuracy order
% forw1: forward or backward
% forw2: forward or backward
% pp1: slope
% pp2: slope
% trace1: input trace
% trace2: input trace
% 
% OUTPUT:
% trace: output trace
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.


n1=size(trace1,1);
nb=2*nw;
eps=e;
eps2=e;

diag=zeros(n1,1);
offd=zeros(n1,nb);



[diag,offd] = regularization(diag,offd,nw,eps,eps2);



[w1,diag,offd] = pwd_define(forw1,diag,offd,n1,nw,pp1);
[w2,diag,offd] = pwd_define(forw2,diag,offd,n1,nw,pp2);

t0=0.5*(trace1(1)+trace2(1));
t1=0.5*(trace1(2)+trace2(2));
t2=0.5*(trace1(n1-1)+trace2(n1-1));
t3=0.5*(trace1(n1)+trace2(n1));

[tmp1] = pwd_set(0,w1,diag,offd,pp1,trace1);
[tmp2] = pwd_set(0,w2,diag,offd,pp2,trace2);

trace=tmp1+tmp2;

trace(1)=trace(1)+eps2*t0;
trace(2)=trace(2)+eps2*t1;
trace(n1-1)=trace(n1-1)+eps2*t2;
trace(n1)=trace(n1)+eps2*t3;

trace = banded_solve(n1,nb,diag,offd,trace);   

return