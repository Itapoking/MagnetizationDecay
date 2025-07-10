function dM = BlochMcConnell(M,  dw, wdip, wrd, R1, R2, M0a, M0b, gridx, gridy, gridz, mask, k)
    % M - vector gridx*gridy*gridz*3, magnetization
    % wdp - vector gridx*gridy*gridz*3, dipolar field
    % wrd - vector gridx*gridy*gridz*3, radiation damping field
    % wrf - vector gridx*gridy*gridz*3, RF field
    % R1, R2, Rrd - constant
    % dwx - vector gridx*gridy*gridz, offset
    dM = ones(size(wdip));
    nstep  = gridx*gridy*gridz;
    Ma = M(1:3*nstep);
    Mb = M(3*nstep + 1:6*nstep);

    dwx = dw(1:nstep) + wrd(1:nstep) + wdip(1:nstep);
    dwy = dw(nstep+1:2*nstep) + wrd(nstep+1:2*nstep) + wdip(nstep+1:2*nstep);
    dwz = dw(2*nstep+1:3*nstep) + wdip(2*nstep+1:3*nstep);% + wrd(2*nstep+1:3*nstep)
    % Ma differential equations: 
    % -k*Ma(1:nstep) + k*Mb(1:nstep) is the exchange ter
    dM(1:nstep) =  dwz.*Ma(nstep+1:2*nstep) - dwy.*Ma(2*nstep+1:3*nstep) - R2*Ma(1:nstep) - k*Ma(1:nstep) + k*Mb(1:nstep);
    dM(nstep + 1:2*nstep) =  dwx.*Ma(2*nstep+1:3*nstep) - dwz.*Ma(1:nstep) - R2*Ma(nstep+1:2*nstep) - k*Ma(nstep+1:2*nstep) + k*Mb(nstep+1:2*nstep);
    dM(2*nstep + 1:3*nstep) =  -dwx.*Ma(nstep+1:2*nstep) + dwy.*Ma(1:nstep) - R1*(Ma(2*nstep+1:3*nstep) - M0a).*mask(:) - k*Ma(2*nstep+1:3*nstep) + k*Mb(2*nstep+1:3*nstep);
    
    % Mb differential equations: 
    dM(3*nstep + 1:4*nstep) =  dwz.*Mb(nstep+1:2*nstep) - dwy.*Mb(2*nstep+1:3*nstep) - R2*Mb(1:nstep) - k*Mb(1:nstep) + k*Ma(1:nstep);
    dM(4*nstep + 1:5*nstep) =  dwx.*Mb(2*nstep+1:3*nstep) - dwz.*Mb(1:nstep) - R2*Mb(nstep+1:2*nstep) - k*Mb(nstep+1:2*nstep) + k*Ma(nstep+1:2*nstep);
    dM(5*nstep + 1:6*nstep) =  -dwx.*Mb(nstep+1:2*nstep) + dwy.*Mb(1:nstep) - R1*(Mb(2*nstep+1:3*nstep) - M0b).*mask(:) - k*Mb(2*nstep+1:3*nstep) + k*Ma(2*nstep+1:3*nstep);
end