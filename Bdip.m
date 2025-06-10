function [Bdipx, Bdipy, Bdipz] = Bdip(M, K, mask, DFTset, gridsize)
    permeability = 1;
    Mx = M(:,:,:,1); My = M(:,:,:,2); Mz = M(:,:,:,3);
    zfx = DFTset(1); zfy = DFTset(1); zfz = DFTset(1);
    gridx = gridsize(1); gridy = gridsize(2); gridz = gridsize(3);

    Mx_padded = padarray(Mx, [zfx, zfy, zfz], 0, 'post');
    My_padded = padarray(My, [zfx, zfy, zfz], 0, 'post');
    Mz_padded = padarray(Mz, [zfx, zfy, zfz], 0, 'post');
    
    Mxk = fftshift(fftn(Mx_padded)); Myk = fftshift(fftn(My_padded)); Mzk = fftshift(fftn(Mz_padded));
    Mxk = (1/6) * (1 - 3*K(:,:,:,1)) .* Mxk * permeability;
    Myk = (1/6) * (1 - 3*K(:,:,:,2)) .* Myk * permeability;
    Mzk = (2/6) * (3*K(:,:,:,3) - 1) .* Mzk * permeability;

    Mxk = ifftshift(Mxk); Myk = ifftshift(Myk); Mzk = ifftshift(Mzk);
    Mx = real(ifftn(Mxk)); My = real(ifftn(Myk)); Mz = real(ifftn(Mzk)); 

    Bdipx = -Mx(1:gridx, 1:gridy, 1:gridz).*mask;
    Bdipy = -My(1:gridx, 1:gridy, 1:gridz).*mask;
    Bdipz = -Mz(1:gridx, 1:gridy, 1:gridz).*mask;

end

