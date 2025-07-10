function [Minit, K, mask, grid, dw] = Initialize_box(gridsize, M0initA, M0initB, shapesize, DFTset, offset, length, height, shim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  gridsize = [nx, ny, nz] number of points for grid 
%  boxsize = [a, b, c] - size of simulation box along Ox, Oy, Oz
%  Minit = [theta, phi, M0, enhancement] - angles for initial magnetization, initial magnetization 
%  shapesize = [L - size of cylinder, R - radius of cylinder]
%  DFTset = [zfx, zfy, zfz] zero filling;
%  center of OxOyOz at the center of grid box
%  shim [z, z2]

    % unpacking values
    gridx = gridsize(1); gridy = gridsize(2); gridz = gridsize(3);
    thetaA = M0initA(1)/180*pi; phiA = M0initA(2)/180*pi; M0A = M0initA(3); enhA = M0initA(4);
    thetaB = M0initB(1)/180*pi; phiB = M0initB(2)/180*pi; M0B = M0initB(3); enhB = M0initB(4);
    zfx = DFTset(1); zfy = DFTset(2); zfz = DFTset(3); 
    L = shapesize(1); R = shapesize(2);
    dx = length/gridx; dy = length/gridy; dz = height/gridz;

    % make offset grid
    x = linspace(-(gridx-1)/2, (gridx-1)/2, gridx);
    y = linspace(-(gridy-1)/2, (gridy-1)/2, gridy);
    z = linspace(-(gridz-1)/2, (gridz-1)/2, gridz);
    grid = zeros(gridx, gridy, gridz, 3);
    [x, y, z] = meshgrid(x, y, z);
    grid(:, :, :, 1) = x; grid(:, :, :, 2) = y; 
    grid(:, :, :, 3) = z;


    dw = ones(gridx, gridy, gridz, 3);
    dw(:,:,:,1) = dw(:,:,:,1)*offset(1);
    dw(:,:,:,2) = dw(:,:,:,2)*offset(2);
    dw(:,:,:,3) = dw(:,:,:,3)*offset(3) + shim(1)/gridz*grid(:,:,:,3) + shim(2)/gridz^2*grid(:,:,:,3).^2;
    dw = dw(:); 
    
    % mask for cylinder
    mask = zeros(gridx, gridy, gridz);
    for k=1:gridz
        if abs(z(1, 1, k)) < L
            for i=1:gridx
                for j=1:gridy
                    r = sqrt(x(i, j, k).^2 + y(i, j, k).^2);
                    if r <= R
                        mask(i, j, k) = 1;
                    end
                end
            end
        end
    end


    % make dimensionless K-space dk = 1/N/dx, dx = 1
    kx = linspace(-1/2, 1/2, gridx + zfx)/dx;
    ky = linspace(-1/2, 1/2, gridy + zfy)/dy;
    kz = linspace(-1/2, 1/2, gridz + zfz)/dz;
    [KX, KY, KZ] = ndgrid(kx, ky, kz);
    Kmag = sqrt(KX.^2 + KY.^2 + KZ.^2);
    KX_unit = KX ./ Kmag;
    KY_unit = KY ./ Kmag;
    KZ_unit = KZ ./ Kmag;
    KX_unit(isnan(KX_unit)) = 0;
    KY_unit(isnan(KY_unit)) = 0;
    KZ_unit(isnan(KZ_unit)) = 0;
    K = zeros(gridx + zfx, gridy + zfy, gridz + zfz, 3);
    K(:,:,:,1) = KX_unit .^ 2;
    K(:,:,:,2) = KY_unit .^ 2;
    K(:,:,:,3) = KZ_unit .^ 2;


    % make M grid array
    MdVA = M0A*enhA;
    Minit = ones(gridx, gridy, gridz, 6);
    Minit(:,:,:,1) = MdVA*sin(thetaA)*cos(phiA)*Minit(:,:,:,1);
    Minit(:,:,:,2) = MdVA*sin(thetaA)*sin(phiA)*Minit(:,:,:,2);
    Minit(:,:,:,3) = MdVA*cos(thetaA)*cos(phiA)*Minit(:,:,:,3);

    Minit(:,:,:,1) = Minit(:,:,:,1).*mask;
    Minit(:,:,:,2) = Minit(:,:,:,2).*mask;
    Minit(:,:,:,3) = Minit(:,:,:,3).*mask;

    MdVB = M0B*enhB;
    Minit(:,:,:,4) = MdVB*sin(thetaB)*cos(phiB)*Minit(:,:,:,1);
    Minit(:,:,:,5) = MdVB*sin(thetaB)*sin(phiB)*Minit(:,:,:,2);
    Minit(:,:,:,6) = MdVB*cos(thetaB)*cos(phiB)*Minit(:,:,:,3);

    Minit(:,:,:,4) = Minit(:,:,:,4).*mask;
    Minit(:,:,:,5) = Minit(:,:,:,5).*mask;
    Minit(:,:,:,6) = Minit(:,:,:,6).*mask;
    
    Minit = Minit(:);
end