% 4/20
% function for HAS that returns wavenumber grids
function [kx2, ky2] = getKgridsquared(Nx,Ny,vox)
    
    % x
    % create wavenumber vector
    if rem(Nx, 2) == 0
        k_vec = ((-Nx/2):(Nx/2-1)) .* 2 * pi ./ (Nx * vox);
    else
        k_vec = (-(Nx-1)/2:(Nx-1)/2) .* 2 * pi ./ (Nx * vox);
    end

    % force middle value to be zero in case 1/Nx is a recurring
    % number and the series doesn't give exactly zero
    k_vec(floor(Nx/2) + 1) = 0;

    % shift wavenumbers to be in the correct order for FFTW
    k_vecx = ifftshift(k_vec);
    
    % y
    % create wavenumber vector
    if rem(Ny, 2) == 0
        k_vec = ((-Ny/2):(Ny/2-1)) .* 2 * pi ./ (Ny * vox);
    else
        k_vec = (-(Ny-1)/2:(Ny-1)/2) .* 2 * pi ./ (Ny * vox);
    end

    % force middle value to be zero in case 1/Nx is a recurring
    % number and the series doesn't give exactly zero
    k_vec(floor(Ny/2) + 1) = 0;

    % shift wavenumbers to be in the correct order for FFTW
    k_vecy = ifftshift(k_vec);

    % create wavenumber grids
    [kx, ky] = meshgrid(k_vecx, k_vecy);
    
    kx2 = kx.^2;
    ky2 = ky.^2;

end