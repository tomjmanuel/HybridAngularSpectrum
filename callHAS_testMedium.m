%% call tom's HAS code
clear all
Nx = 128; Ny = 128; Nz = 128;

% create Amps and phases
zp = 30; %zero padding
Amps = ones([Nx Ny]);
Amps(1:zp,:)=0;
Amps(end-zp:end,:)=0;
Amps(:,1:zp)=0;
Amps(:,end-zp:end)=0;

phases = zeros([Nx Ny]);

% frequency
f = 1.1E6; %1MHz

% voxel size (isotropic for now)
vox = 0.25E-3; %m

% medium properties (homogenous for now)
sos = 1480.*ones([Nx Ny Nz]);
density = 997.*ones([Nx Ny Nz]);

load('chunk.mat')
sos(bar)=3000;
density(bar)=1500;
%sos(bar)=3000;
%density(bar)=1100;

p0 = complex(Amps,atan(phases));

figure
%additional zero padding
zpf=64;
p0 = padarray(p0,[zpf zpf]);
sos = padarray(sos,[zpf zpf 0],1480);
density = padarray(density,[zpf zpf 0],997);

[p, pr, pf]  = HASv4(p0,vox,f,sos,density);

% remove padding
p = p(zpf+1:end-zpf,zpf+1:end-zpf,:);

%p = HASv2(p0,vox,f,sos,density);
%subplot(131)
imagesc(squeeze(abs(p(:,64,:))))
axis image
title('p total')
% subplot(132)
% imagesc(squeeze(abs(pr(:,64,:))))
% title('p ref')
% axis image
%  subplot(133)
% imagesc(squeeze(abs(pf(:,64,:))))
% title('p forward')   
% axis image
