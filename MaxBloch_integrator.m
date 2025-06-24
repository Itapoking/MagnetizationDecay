function dM = MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0, gridx, gridy, gridz, mask)
    % M - vector gridx*gridy*gridz*3, magnetization
    % wdp - vector gridx*gridy*gridz*3, dipolar field
    % wrd - vector gridx*gridy*gridz*3, radiation damping field
    % wrf - vector gridx*gridy*gridz*3, RF field
    % R1, R2, Rrd - constant
    % dwx - vector gridx*gridy*gridz, offset
    dM = ones(size(wdip));
    nstep  = gridx*gridy*gridz;

    dwx = dw(1:nstep) + wrd(1:nstep) + wdip(1:nstep);
    dwy = dw(nstep+1:2*nstep) + wrd(nstep+1:2*nstep) + wdip(nstep+1:2*nstep);
    dwz = dw(2*nstep+1:3*nstep) + wdip(2*nstep+1:3*nstep);% + wrd(2*nstep+1:3*nstep)
    
    dM(1:nstep) =  dwz.*M(nstep+1:2*nstep) - dwy.*M(2*nstep+1:3*nstep) - R2*M(1:nstep);
    dM(nstep+1:2*nstep) =  dwx.*M(2*nstep+1:3*nstep) - dwz.*M(1:nstep) - R2*M(nstep+1:2*nstep);
    dM(2*nstep+1:3*nstep) =  -dwx.*M(nstep+1:2*nstep) + dwy.*M(1:nstep) - R1*(M(2*nstep+1:3*nstep) - M0).*mask(:);

end