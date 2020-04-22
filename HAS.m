%HAS v4 has reflected pressures. Pretty sure indexing is all correct.
% wondering if zero padding more would help. 

% HASv3 uses weighted average method developed by steven leung from
% stanford
% HAS v3 doesn't have reflected pressures
% HAS v3 may not have indexing in z correct
% base function of HAS is to propogate a plane forward in space
% inputs: 
%   amps    [mxn]
%   phases  [mxn] [radians]
%   freq cyc/s 
%   vox:    dx (isotropic for now)
%   sos     [mxnxp] 
%   density [mxnxp]
%   alpha   %attenuation, leave out for now

% outputs:
%   p: total pressure
%   Refp: reflected pressure
%   Forp: forward pressure

function [p, Refp, Forp] = HAS(p0,vox,f,sos,density)
    
    dim = size(sos);
    Nx = dim(1);
    Ny = dim(2);
    Nz = dim(3);
    dz = vox;

    % allocate pfield memory
    p = zeros([Nx Ny Nz]);
    p = complex(p,0);
    p(:,:,1) = p0;
    
    % create wavenumber grids
    [kx2,ky2] = getKgridsquared(Nx,Ny,vox);

    % calculate reflection and transmission coefficients
    imped = density.*sos; % impedance
    Ref = (imped(:,:,2:end)...
        -imped(:,:,1:end-1))./(imped(:,:,2:end)+imped(:,:,1:end-1)); %reflection coeff
    T = 1+Ref; % transmission coeff
    %Ref(z==1) refers to reflection coefficient btwn layers 2&1
    
    % allocate memory for reflected propogations
    Refp = p; %reflected pressure
    
    % propogate forward
    for z=1:Nz-1

        % get wavenumber matrix for given layer
        kLayer = 2*pi*f./sos(:,:,z);
        
        % get kprime forward
        kpf = weightedAverageSteve(abs(p(:,:,z)), kLayer);
        
        % cast pressure into f domain (eqn 8)
        ppf = fft2(p(:,:,z));
        
        % propogate in spatial frequency domain (eqn 9)
        kz = sqrt(kpf.^2 - kx2 - ky2);
        H = exp(1i.*dz.*kz);
        ppf = ifft2(ppf.*H);
        
        % propogate spatial domain (eqn 7)
        % T scales amps based on impedance mismatch
        rp = dz; % Steve also didn't get the r prime bit
        dbn = kLayer - kpf; % difference in spatial properties of medium
        p(:,:,z+1) = T(:,:,z).*ppf.*exp(rp*1i.*dbn); %pn'
        
        % store refelcted pressure at this plane (eqn 11)
        Refp(:,:,z) = Ref(:,:,z).*p(:,:,z+1); %Ref(z) is actually Ref(z+1) due to line 38 syntax
    end
    Forp = p;
    
    % backwards propogate and add in reflections at each layer
    % in prop forward, when z=1, calc for z=2
    % here, z=2, calc for 1
    % Refp(:,:,n) gives reflection that happened from plane n-1 to n
    % inject them at n-1?
    for z=Nz:-1:2
        % skip this layer if there is nothign in refP here
        if ~(sum(sum(Refp(:,:,z)))==0)
            % get wavenumber matrix for given layer
            kLayer = 2*pi*f./sos(:,:,z);

            % get kprime forward
            kpf = weightedAverageSteve(abs(Refp(:,:,z)), kLayer);

            % cast pressure into f domain (eqn 8)
            ppf = fft2(Refp(:,:,z));

            % propogate in spatial frequency domain (eqn 9)
            kz = sqrt(kpf.^2 - kx2 - ky2);
            H = exp(1i.*dz.*kz);
            ppf = ifft2(ppf.*H);

            % propogate spatial domain (eqn 7)
            % T scales amps based on impedance mismatch
            rp = dz; % Steve also didn't get the r prime bit
            dbn = kLayer - kpf; % difference in spatial properties of medium
            Refp(:,:,z-1) = Refp(:,:,z-1)+ppf.*exp(rp*1i.*dbn); %pn'

        end
        
    end
    
    %sum transmitted and reflected fields 
    p = Forp + Refp;
end