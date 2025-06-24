function wdip = Bfield_DP(M, K, mask, gridsize, DFTset, gamma, mu0, length, height)
    % some explanation
    % Kernel operator is dimensionless
    % wdf  = gamma*Bdf = gamma*mu0*dMdf, where dMdf is scaled magnetization
    % of primitive cell of the grid. It is already scaled by initial
    % construction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%|  INPUT   |%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  gridsize = [nx, ny, nz] number of points for grid 
    %  M = Mx, My, Mz - magnetization at each grid point
    %  K = Kx, Ky, Kz - inverse space vector
    %  shapesize = [L - size of cylinder, R - radius of cylinder]
    %  DFTset = [zfx, zfy, zfz] zero filling;
    %  center of OxOyOz at the center of grid box
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = reshape(M, [size(mask) 3]);
    Mx = M(:,:,:,1); My = M(:,:,:,2); Mz = M(:,:,:,3);
    zfx = DFTset(1); zfy = DFTset(2); zfz = DFTset(3);
    gridx = gridsize(1); gridy = gridsize(2); gridz = gridsize(3);
    Vgrid = gridx*gridy*gridz;
    wdip = zeros(gridx, gridy, gridz, 3);
     
    Mxaver = (1/3)*mean(Mx, "All")*Vgrid/sum(mask(:,:,:),"All");
    Myaver = (1/3)*mean(My, "All")*Vgrid/sum(mask(:,:,:),"All");
    Mzaver = (-2/3)*mean(Mz, "All")*Vgrid/sum(mask(:,:,:),"All");
    
    Mx_padded = padarray(Mx, [zfx, zfy, zfz], 0, 'post');
    My_padded = padarray(My, [zfx, zfy, zfz], 0, 'post');
    Mz_padded = padarray(Mz, [zfx, zfy, zfz], 0, 'post');
    
    Mxk = fftshift(fftn(Mx_padded)); Myk = fftshift(fftn(My_padded)); Mzk = fftshift(fftn(Mz_padded));
    Mxk = (1/6) * (1 - 3*K(:,:,:,1)) .* Mxk;
    Myk = (1/6) * (1 - 3*K(:,:,:,2)) .* Myk;
    Mzk = (2/6) * (3*K(:,:,:,3) - 1) .* Mzk;

    Mxk = ifftshift(Mxk); Myk = ifftshift(Myk); Mzk = ifftshift(Mzk);
    Mx = real(ifftn(Mxk)); My = real(ifftn(Myk)); Mz = real(ifftn(Mzk)); 

    wdip(:, :, :, 1) = (Mxaver - Mx(1:gridx, 1:gridy, 1:gridz)).*mask;
    wdip(:, :, :, 2) = (Myaver - My(1:gridx, 1:gridy, 1:gridz)).*mask;
    wdip(:, :, :, 3) = (Mzaver - Mz(1:gridx, 1:gridy, 1:gridz)).*mask;

    wdip(:, :, :, 1) = -Mx(1:gridx, 1:gridy, 1:gridz).*mask;
    wdip(:, :, :, 2) = -My(1:gridx, 1:gridy, 1:gridz).*mask;
    wdip(:, :, :, 3) = -Mz(1:gridx, 1:gridy, 1:gridz).*mask;

    wdip = 2*pi*mu0*gamma*wdip(:);
end

 