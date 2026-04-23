
% Example for running RTM.m

clear all;clc; close all;

% input data, t, tmax
load('data.mat')
data = Ydata;
data = BGRm(data);
t = t.*1e-9;
tmax = t(1024).*1e-9;

% cut delay time
n_c = 152;
data(1:n_c,:) = []; 

%% set parameter of model 
% set x, z, dx, dz
x = 1.8; 
z = 1.0;
dx = 0.005;
dz = 0.005;
x = 0:dx:x;
z = 0:dz:z;
nx = size(x,2);
nz = size(z,2);

% set ep(relative permittivity), mu(permeability), sig(conductivity)
ep = zeros(nz,nx); ep(1:10,:) = 1; ep(11:end,:) = 2.4;
mu = zeros(nz,nx); mu(:) = 1;
sig = zeros(nz,nx); sig(:) = 0.001;

% imagesc(ep)

ep = 4.*ep; % 1/2 velocity
ep = ep';
mu = mu';
sig = sig';

% interpolate electrical property grids to proper spatial discretization

x2 = min(x):dx/2:max(x);
z2 = min(z):dx/2:max(z);
ep2 = gridinterp(ep,x,z,x2,z2,'nearest');
mu2 = gridinterp(mu,x,z,x2,z2,'nearest');
sig2 = gridinterp(sig,x,z,x2,z2,'nearest');
%% set time

dt= (t(2)-t(1));  %%% time step    
t(:,length(t)-n_c:end) = [];

% pad electrical property matrices for PML absorbing boundaries
npml = 40;  % number of PML boundary layers
[ep3,x3,z3] = padgrid(ep2,x2,z2,2*npml);
[mu3,x3,z3] = padgrid(mu2,x2,z2,2*npml);
[sig3,x3,z3] = padgrid(sig2,x2,z2,2*npml);

%% create source and receiver location matrices
srcx = (0:0.025:1.8)';
srcz = 0.0*ones(size(srcx));
recx = srcx + 0.000;
recz = srcz;
srcloc = [srcx srcz];
recloc = [recx recz];

% run
tic;
[R_pro,srcx,srcz,recx,recz] = GPR_RTM(ep3,mu3,sig3,x3,z3,srcloc,recloc,t,npml,data);
disp(' ');
disp(['Total running time = ',num2str(toc/3600),' hours']);
[nz,nx] = size(R_pro);

figure;
imagesc(x,z,data);
colormap('gray');
ylabel('Depth (m)'); 
xlabel('Distance (m)');
set(gca,'FontSize',14);
colorbar
pbaspect([2 1 1])
set(gca,'ytick',0:0.2:1)

figure; 
imagesc(x,z,R_pro(npml:nz-npml,npml:nx-npml));
% caxis([-1e-1 1e-1]);
colormap 'gray';
ylabel('Depth (m)'); 
xlabel('Distance (m)');
set(gca,'FontSize',14);
colorbar
pbaspect([2 1 1])
set(gca,'ytick',0:0.2:1)


