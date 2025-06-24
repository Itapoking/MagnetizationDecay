function wrd = Bfield_RD(M, psi, mu0, gamma, nuQ, mask)
    % psi - angle between RD field and magnetization
    % nuQ - coil constant
    % muo, gamma - vacuu magnetic permeability and gyromagnetic ration

    [Mx, My, Mz, Mnorm] = M_averaging(M, mask);
    wrd = ones([size(mask) 3]);
    psi = psi/180*pi;
    sizemask = size(mask);
    Bx_RD = 1*(-Mx*sin(psi) - My*cos(psi));
    By_RD = -1*(-Mx*cos(psi) + My*sin(psi));
    Bz_RD = 0;


    wrd(:,:,:,1) = Bx_RD*wrd(:,:,:,1).*mask;
    wrd(:,:,:,2) = By_RD*wrd(:,:,:,2).*mask;
    wrd(:,:,:,3) = Bz_RD*wrd(:,:,:,3).*mask;
    wrd = 2*pi*nuQ*mu0*gamma/2*wrd(:);
end
% we
% A = sin
% B = cos
% sin Mx - cos My
% cos Mx + sin My