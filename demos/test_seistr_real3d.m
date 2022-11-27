%% Demo for structure-oriented filtering algorithm of 3D data
%
%  The input free parameters of the slope estimation and structural
%  filtering are listed as follows: 
%-------------------------------------------------------------------------
%  1. slope estimation (str_dip2d.m)
%  din: input data for slope estimation
%  niter: number of non-linear iterations (default value 5-10)
%  liter: number of linear iterations (default value 10-20)
%  order: accuracy order of PWD filter (default value 1 or 2)
%  eps_dv: regualrization parameter in the non-linear iteration (default value 0.01)
%  eps_cg: regualrization parameter in the linear iteration (CG) (default value 1)
%  eps_cg: tolerance in the linear iteration (CG) (default value 0.000001)
%  rect:  size of the triangle smoothing operator (default value 5-10)
%  verb: verbosity flag that controls if printing the iteration number of CG (default value 1)
%-------------------------------------------------------------------------
%  2. structural filtering (str_pwsmooth_lop2d.m)
%  dn: noisy data
%  dipi: inline slope field
%  dipx: inline slope field
%  r1,r2: spray radius (smoothing length) (default value 1-4)
%  order: accuracy order of PWD filter (default value 1 or 2)
%  eps: regularization parameter (default value 0.01);

%% load 3D data
clear;clc;close;
addpath(genpath('../seistr'))

%The input 3D source data file "real3d.bin" can be downloaded from
%https://github.com/chenyk1990/reproducible_research/blob/master/drr3d/matfun/real3d.bin,

fid=fopen('data/real3d.bin','r');
d=fread(fid,[300,1000],'float');
d=reshape(d,300,100,10);
d=d(200:300,50:100,:);
cmp=d./max(max(max(d)));

cmpn=cmp;
%% 3D slope calculation (inline and xline)
% default parameter values are suitable for most cases
[dipi,dipx] = str_dip3d(cmpn);

%% Structural smoothing
r1=2;
r2=2;
eps=0.01;
order=2;

cmpn_d1=str_pwsmooth_lop3d(cmpn,dipi,dipx,r1,r2,eps,order);

%% conventional
temp=cmpn*0;%temp is the smoothed data by the conventional method
for i=1:size(cmpn,1)
    for j=1:size(cmpn,2)
        temp(i,j,:)=smooth(cmpn(i,j,:));
    end
end

for i=1:size(cmpn,1)
    for j=1:size(cmpn,3)
        temp(i,:,j)=smooth(temp(i,:,j));
    end
end

%% visualization
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(3,2,1);
str_plot3d(cmpn,[50,20,5],[1:101],[1:51],[1:10]);caxis([-0.5,0.5]);zlim([1,101]);xlim([1,51]);ylim([1,10]);title('Noisy');
subplot(3,2,3);
str_plot3d(temp,[50,20,5],[1:101],[1:51],[1:10]);caxis([-0.5,0.5]);zlim([1,101]);xlim([1,51]);ylim([1,10]);title('Denoised (Mean)');
subplot(3,2,4);
str_plot3d(cmpn-temp,[50,20,5],[1:101],[1:51],[1:10]);caxis([-0.5,0.5]);zlim([1,101]);xlim([1,51]);ylim([1,10]);title('Noise (Mean)');
subplot(3,2,5);
str_plot3d(cmpn_d1,[50,20,5],[1:101],[1:51],[1:10]);caxis([-0.5,0.5]);zlim([1,101]);xlim([1,51]);ylim([1,10]);title('Denoised (SOMEAN)');
subplot(3,2,6);
str_plot3d(cmpn-cmpn_d1,[50,20,5],[1:101],[1:51],[1:10]);caxis([-0.5,0.5]);zlim([1,101]);xlim([1,51]);ylim([1,10]);title('Noise (SOMEAN)');
print(gcf,'-dpng','-r300','test_seistr_real3d.png');

