function [Minit, K, mask,x,y,z] = Mtest(gridsize, boxsize, Minit, shapesize, DFTset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  gridsize = [nx, ny, nz] number of points for grid 
%  boxsize = [a, b, c] - size of simulation box along Ox, Oy, Oz
%  Minit = [theta, phi, M0] - angles for initial magnetization, initial magnetization 
%  shapesize = [L - size of cylinder, R - radius of cylinder]
%  DFTset = [zfx, zfy, zfz] zero filling;
%  center of OxOyOz at the center of grid box
%%%%%%%%%%% Example %%%%%%%%%%%%%% 
 % gridsize = [16, 16, 16];
 % boxsize = [2, 2, 10];
 % Minit = [90, 0, 1];
 % shapesize = [5, 0.5];
 % DFTset = [100, 100, 100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gridx = gridsize(1); gridy = gridsize(2); gridz = gridsize(3);
    a = boxsize(1); b = boxsize(2); c = boxsize(3);
    theta = Minit(1)/180*pi; phi = Minit(2)/180*pi; M0 = Minit(3);
    zfx = DFTset(1); zfy = DFTset(1); zfz = DFTset(1); 
    L = shapesize(1); R = shapesize(2);

    % make grid box it is dimensionless
    % step length is 1, since number of points is N
    x = linspace(-(gridx-1)/2, (gridx-1)/2, gridx); %da = a/(gridx); x = da*x;
    y = linspace(-(gridy-1)/2, (gridy-1)/2, gridy); %db = b/(gridy); y = db*y;
    z = linspace(-(gridz-1)/2, (gridz-1)/2, gridz); %dc = c/(gridz); z = dc*z;
    [x, y, z] = meshgrid(x, y, z);


    % make dimensionless K-space dk = 1/N/dx, dx = 1
    kx = linspace(-1/2, 1/2, gridx + zfx);
    ky = linspace(-1/2, 1/2, gridy + zfy);
    kz = linspace(-1/2, 1/2, gridz + zfz);   

    % Трёхмерная решётка
    [KX, KY, KZ] = ndgrid(kx, ky, kz);  % аналог indexing='ij' в NumPy
    
    % Модуль вектора
    Kmag = sqrt(KX.^2 + KY.^2 + KZ.^2);
    
    % Нормированные компоненты, обработка деления на ноль
    KX_unit = KX ./ Kmag;
    KY_unit = KY ./ Kmag;
    KZ_unit = KZ ./ Kmag;
    
    KX_unit(isnan(KX_unit)) = 0;
    KY_unit(isnan(KY_unit)) = 0;
    KZ_unit(isnan(KZ_unit)) = 0;
    
    % Создание массива K размером (nx+paddx, ny+paddy, nz+paddz, 3)
    K = zeros(gridx + zfx, gridy + zfy, gridz + zfz, 3);
    
    K(:,:,:,1) = KX_unit .^ 2;
    K(:,:,:,2) = KY_unit .^ 2;
    K(:,:,:,3) = KZ_unit .^ 2;


    % mask for cylinder
    mask = zeros(gridx, gridy, gridz);
    for k=1:gridz
        for i=1:gridx
            for j=1:gridy
                r = sqrt(x(i, j, 1).^2 + y(i, j, 1).^2);
                if r <= gridx/3
                    mask(i, j, k) = 1;
                end
            end
        end
    end
    
    % make M grid array
    Minit = M0*ones(gridx, gridy, gridz, 3);
    Minit(:,:,:,1) = M0*sin(theta)*cos(phi)*Minit(:,:,:,1);
    Minit(:,:,:,2) = M0*sin(theta)*sin(phi)*Minit(:,:,:,2);
    Minit(:,:,:,3) = M0*cos(theta)*cos(phi)*Minit(:,:,:,3);

    Minit(:,:,:,1) = Minit(:,:,:,1).*mask;
    Minit(:,:,:,2) = Minit(:,:,:,2).*mask;
    Minit(:,:,:,3) = Minit(:,:,:,3).*mask;
end