function [Minit, K, mask, grid, dw] = Initialize_box(gridsize, M0init, shapesize, DFTset, offset, length, height, shim)
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
    theta = M0init(1)/180*pi; phi = M0init(2)/180*pi; M0 = M0init(3); enh = M0init(4);
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
    MdV = M0*enh;
    Minit = ones(gridx, gridy, gridz, 3);
    Minit(:,:,:,1) = MdV*sin(theta)*cos(phi)*Minit(:,:,:,1);
    Minit(:,:,:,2) = MdV*sin(theta)*sin(phi)*Minit(:,:,:,2);
    Minit(:,:,:,3) = MdV*cos(theta)*cos(phi)*Minit(:,:,:,3);

    Minit(:,:,:,1) = Minit(:,:,:,1).*mask;
    Minit(:,:,:,2) = Minit(:,:,:,2).*mask;
    Minit(:,:,:,3) = Minit(:,:,:,3).*mask;
    Minit = Minit(:);
end