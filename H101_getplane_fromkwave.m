% Tom
% Simulate idea phase of H101

%% Setup Sim space Extent
clear all
Nx = 256; Ny = 256; Nz = 256;

% using h115, focusing down x dim
XrangeSim = .08;
YrangeSim = .08;
ZrangeSim = .08;

dxSim = XrangeSim/Nx;
dySim = YrangeSim/Ny;
dzSim = ZrangeSim/Nz;

%% setup transducer
% H115
transducer_radius_m = (63.2e-03);
transducer_height_m = (63.2 - 54.5)*1e-03;
transducer_radius_pix = round( transducer_radius_m / dySim );
transducer_height_pix = round( transducer_height_m / dxSim );
sphere = makeSphericalSection( transducer_radius_pix, transducer_height_pix );

% trasnducer has different size but same dx/dy/dz as simspace
dimTrans = size(sphere);

%% Place transducer at top of simulation space

transGrid = zeros([Nx Ny Nz]);
xoffset = 6; % pix
yoffset = 15;
zoffset = 15;

transGrid(xoffset:dimTrans(1)+xoffset-1,...
    yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1) = sphere;

% store focus location
xF = xoffset + transducer_radius_pix;
yF = yoffset + floor(dimTrans(2)/2);
zF = zoffset + floor(dimTrans(3)/2);

% have plane be at flat face of xdcr (will be lowest pixel after rot)
xLoc = round(xoffset + transducer_height_m/dxSim);
fLoc = [xF yF zF];

%% make grid that will be converted to nifti (makeNiftiTransducer.m)
xdcr = transGrid.*255;
bs = 3; % add points around focus to make vis easier
xdcr(fLoc(1)-bs:fLoc(1)+bs,fLoc(2)-bs:fLoc(2)+bs,fLoc(2)-bs:fLoc(2)+bs)=255;
xdcr(fLoc(1),fLoc(2),fLoc(3))=1; % set focus to one for future ref
xdcr(xoffset,fLoc(2),fLoc(3))=2; % set one pixel at top center = 2 for future ref

%% crop xdcr to match extent
xdcr = xdcr(xoffset:xoffset+xF+bs,yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1);
H101mask = xdcr;
save('H101mask.mat','H101mask');

%% setup medium
medium.sound_speed = ones([Nx Ny Nz])*1480;
medium.density = ones([Nx Ny Nz])*997;

% Finally create the full grid
kgrid = kWaveGrid(Nx, dxSim, Ny, dySim, Nz, dzSim);
[kgrid.t_array, dt_SIM] = makeTime(kgrid, medium.sound_speed);
Nt = length(kgrid.t_array);

%% setup sensor
% record orthoganal slices at focus
sensor.mask = zeros([Nx Ny Nz]);
sensor.mask(xLoc,:,:)=1;
sensor.record={'p'};

%% setup source
source.p_mask = transGrid;
% create source vector 
%freq = 1.1*1E6; %Mhz
freq = 1.1*1E6; 

% Aubrey used 2 cycle sinusoid burst with guass apodization
Ncyc = 20;
t = kgrid.t_array;
dt = kgrid.dt;
nptspulse= round( Ncyc * 1 / (dt*freq) ); % npts in our pulse
tp = dt.*(1:1:nptspulse); % time vector for pulse
pulse = sin(2*pi*freq*tp); %pure ncycle sin wave

win = gausswin(nptspulse); %guassian window
pvec = zeros([length(t) 1]);

amp = 1;
pvec(1:nptspulse)=amp.*pulse.*win';
source.p = pvec';

%% run the thing
%input_args = {'PlotSim', true, 'DisplayMask', source.p_mask, 'DataCast',
%'gpuArray-single', 'PMLInside', false, 'PlotPML', false};
comptype='gpuArray-single';

%run the simulation
[sensor_data] = kspaceFirstOrder3DG(kgrid, medium, source, sensor);

%% reshape sensor_data
pout = zeros([Ny Nz Nt]);
curr=1;
for k=1:Nz
    for j=1:Ny
        for i=1:Nx
            if sensor.mask(i,j,k)==1
                pout(j,k,:)=sensor_data.p(curr,:);
                curr=curr+1;
            end
        end
    end
end

