clc
clear all

%  tsoft/Bsoft and everything from shapereader() it is selective pulse
%  which Rob aws using. time step dt in diff equations is fixed because it
%  comes form selective pulse

% shim should be [0 0] for uniform magnetic field in the spectrometer

% to change RD field you can change parameter nuQ (which is from 1 to 100)
% to change DP field you can change parameter enhancdip (it is just scale
% constant for dipolar field which is calculated in Bfield_DP
% M0 vector it is initial conditions for your simulation. theta and phi
% define magnetization in each cell. M0(3) is absolute value of water
% magnetization at our condition. M0(4) is enhancement factor which we have
% in dissolution (from 10 to 50, but more probably around 20)


% simulation parameters. Long cylinder
% you can increase number of points in the mesh (pointx and pointz) if it
% changes your result (you just make more precise mesh)
length = 2.5/1000; height = 10*2.5/1000; % height and length of simulation box in mm
pointx = 20; pointz = 20;
gridsize = 1*[pointx, pointx, pointz];
DFTset = 2*gridsize;                % Zero padding for FFT. 2x of mesh dimension should be fine
shapesize = [pointz/1.1, pointx/3]; % which part of mesh use for cylinder

shim = [0*2*pi 0*2*pi]; % Z and Z2 gradients. For uniform spectrometer field should be 0 and 0
tsoft = 3237*1e-6; % total duration of selective pulse
Bsoft = 1896*2*pi; % Max B1 field of selective pulse
[tpul, Bampx, Bampy, points, amplitude, phase] = shapereader("./reburp1000.xlsx", tsoft, Bsoft);  % 1000 points tpul - time vector of selective pulse; Bampx - vector of B1 amplitude of selective pulse

% magnetization parametes
% |M0|in A/m units; enhancement = 1 for thermal and around 20 for
% hyperpolarized sample
M0 = [0, 0, 0.038, 20];  % theta, phi, |M0|, enhancement
offset = 1*[0, 0, -1776*2*pi]; % -1776*2*pi it is frequency difference between water signal and pulse frequency (you can change +/- sign)
R1 = 1/12; R2 = 1/2; % R1 and R2 is it 1/T1 and 1/T2 which is fixed and you don't need to change it.


% enhancdip it is scale factor of dipolar field. I'm still not sure about
% real values but looks like 
% enhancdip = 1 gives normal values around 4 Hz for thermal M0.
% I think you can change it from 1 to 60
enhancdip=40;

% circuit parameters nuQ it is parameter which define strenght of RD field.
% Normal values of nuQ is 50. But you can change from 20 to 60.
% psi if is "Impedance shift" which is not more than 1 degree in our case.
% You can put 0 and forget about it.
nuQ = 43.57*0.5; psi = 0;

% constant
mu0 = 4*pi*1e-7; gamma = 42.57*1e6;

[Minit, K, mask, grid, dw] = Initialize_box(gridsize, M0, shapesize, DFTset, offset, length, height, shim);
dt = tpul(1);
Nstep = 1e3;
observables = zeros([12*Nstep 5]);
fields = zeros([12*Nstep 7]);
nstepsmal = gridsize(1)*gridsize(2)*gridsize(3); 
% wdip = 0*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0, length, height);
% wrd  = 0*Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);

% wrd1 = reshape(wrd, [size(mask) 3]);
% plot(wrd1(:,:,pointz/2,2))
%  
% % 
for n = 1:Nstep
    n
    dw(1:nstepsmal) = Bampx(n)*ones([nstepsmal 1]); % here I put B1 values to x field from selective pulse
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    [wdipx, wdipy, wdipz, wdipnorm] = M_averaging(wdip, mask);
    [wrdx, wrdy, wrdz, wrdnorm] = M_averaging(wrd, mask);
    observables(n, 1) = dt*(n-1);
    observables(n, 2) = Mx; observables(n, 3) = My; observables(n, 4) = Mz; observables(n, 5) = Mnorm;
end
% % 
% % % % dt = 0.0001
% % % % % dt=2*dt;

% % dt = 2*dt;
for n = 1:2*Nstep
    n + 1000
    dw(1:nstepsmal) = zeros(nstepsmal,1);%Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    observables(n + 1000, 1) = dt*(n + 1000-1);
    observables(n + 1000, 2) = Mx; observables(n + 1000, 3) = My; observables(n + 1000, 4) = Mz; observables(n + 1000, 5) = Mnorm;
end

for n = 1:Nstep
    n + 3000
    dw(1:nstepsmal) = Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    [wdipx, wdipy, wdipz, wdipnorm] = M_averaging(wdip, mask);
    [wrdx, wrdy, wrdz, wrdnorm] = M_averaging(wrd, mask);
    observables(n+3000, 1) = dt*(n+3000-1);
    observables(n+3000, 2) = Mx; observables(n+3000, 3) = My; observables(n+3000, 4) = Mz; observables(n+3000, 5) = Mnorm;
    % fields(n, 1) = dt*(n-1);
    % fields(n, 2) = wdipx; fields(n, 3) = wdipy; fields(n, 4) = wdipz;
    % fields(n, 5) = wrdx; fields(n, 6) = wrdy; fields(n, 7) = wrdz;
end

for n = 1:2*Nstep
    n + 4000
    dw(1:nstepsmal) = zeros(nstepsmal,1);%Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    observables(n+4000, 1) = dt*(n+4000-1);
    observables(n+4000, 2) = Mx; observables(n+4000, 3) = My; observables(n+4000, 4) = Mz; observables(n+4000, 5) = Mnorm;
end

for n = 1:Nstep
    x = n + 6000
    dw(1:nstepsmal) = Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    [wdipx, wdipy, wdipz, wdipnorm] = M_averaging(wdip, mask);
    [wrdx, wrdy, wrdz, wrdnorm] = M_averaging(wrd, mask);
    observables(x, 1) = dt*(x-1);
    observables(x, 2) = Mx; observables(x, 3) = My; observables(x, 4) = Mz; observables(x, 5) = Mnorm;
end

for n = 1:2*Nstep
    x = n + 7000
    dw(1:nstepsmal) = zeros(nstepsmal,1);%Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    observables(x, 1) = dt*(x-1);
    observables(x, 2) = Mx; observables(x, 3) = My; observables(x, 4) = Mz; observables(x, 5) = Mnorm;
end
for n = 1:Nstep
    x = n + 9000
    dw(1:nstepsmal) = Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    [wdipx, wdipy, wdipz, wdipnorm] = M_averaging(wdip, mask);
    [wrdx, wrdy, wrdz, wrdnorm] = M_averaging(wrd, mask);
    observables(x, 1) = dt*(x-1);
    observables(x, 2) = Mx; observables(x, 3) = My; observables(x, 4) = Mz; observables(x, 5) = Mnorm;
end

for n = 1:2*Nstep
    x = n + 10000
    dw(1:nstepsmal) = zeros(nstepsmal,1);%Bampx(n)*ones([nstepsmal 1]);
    wdip = enhancdip*Bfield_DP(Minit, K, mask, gridsize, DFTset, gamma, mu0);
    wrd  = Bfield_RD(Minit, psi, mu0, gamma, nuQ, mask);
    [t, M] = ode113(@(t, M) MaxBloch_integrator(M, dw, wdip, wrd, R1, R2, M0(end-1), gridsize(1), gridsize(2), gridsize(3), mask), linspace(0, dt, 3), Minit);
    Minit = M(3, :);
    [Mx, My, Mz, Mnorm] = M_averaging(Minit, mask);
    observables(x, 1) = dt*(x-1);
    observables(x, 2) = Mx; observables(x, 3) = My; observables(x, 4) = Mz; observables(x, 5) = Mnorm;
end


% % % 
% % % N = size(observables(:, 2),1);
% % % vNyq = 1/2/dt;
% % % dv = (N-1)/(N*observables(end, 1));
% % % omega = -vNyq + linspace(0, N-1, N)*dv;
% % % plot(omega/500, real(fftshift(fft(observables(:, 2)+1i*observables(:, 3)))))
% % % xlim([-20 1])
% % % 
% 1-3 4-6 7-9 10-12
ncut=4000;
sig = observables(ncut:(ncut+2000),2) + 1i*observables(ncut:(ncut+2000),3);
N  = size(sig,1);  
vNyq = 1/(2*dt);      
zfFactor = 8;            
N_zf = zfFactor * N;        
sig = [sig; zeros(N_zf-N,1)];
dv = 1/(N_zf*dt);          
omega = -vNyq + (0:N_zf-1).' * dv;
spectrum = fftshift(fft(sig));
plot(omega/500, real(spectrum));
xlim([-20 20])

% % % 
% M12 = reshape(Minit, size(grid));
% 
% plot(squeeze(grid(20,20,:,3))*height/40*100,squeeze(M12(pointx/2,pointx/2,:,3)),'bo-','LineWidth',2)
% xlabel('Z axis, cm','FontSize',14);
% ylabel('W_d_i_p, Hz','FontSize',14);
% title('DDF field strength, H/d = 10','FontSize',14);
% set('LineWidth',4, 'FontName','Helvetica','FontSize',11);

% plot(squeeze(M12(4,4,:,3)))
 plot(observables(:,1), observables(:,2),'bo-',observables(:,1), observables(:,3),'ro-',observables(:,1), observables(:,4),'o-',observables(:,1), observables(:,5),'b-')
% 
% subplot()
% plot(fields(:,1), fields(:,2),'bo-',fields(:,1), fields(:,3),'ro-',fields(:,1), fields(:,4),'ko-')
% plot(fields(:,1), fields(:,5),'bo-',fields(:,1), fields(:,6),'ro-',fields(:,1), fields(:,7),'ko-')
% 
% % 
% M = reshape(Minit, size(grid));
% nplx=1;
% nplz=1;
% quiver3(grid(1:nplx:end,1:nplx:end,1:nplz:end,1)*length/20*100, grid(1:nplx:end,1:nplx:end,1:nplz:end,2)*length/20*100, grid(1:nplx:end,1:nplx:end,1:nplz:end,3)*height/40*100, M12(1:nplx:end,1:nplx:end,1:nplz:end,1), M12(1:nplx:end,1:nplx:end,1:nplz:end,2), M12(1:nplx:end,1:nplx:end,1:nplz:end,3),'AutoScale','on','AutoScaleFactor',0.8, 'LineWidth',1.1,'Color',[0 0.447 0.741]);              % MATLAB blue
% box on; 
% % axis equal tight;
% xlabel('X axis, cm','FontSize',12);                     % подписи
% ylabel('Y axis, cm','FontSize',12);
% zlabel('Z axis, cm','FontSize',12);
% title('3‑D DDF field strength, Hz','FontSize',14);
% set(gca, 'LineWidth',1, 'FontName','Helvetica','FontSize',11,'TickDir','out');
% ax = gca;
% ax.DataAspectRatio=[1 1 4]
% % % plot(omega(1001:2000), -real(fftshift(fft(observables(1001:2000, 2)+1i*observables(1001:2000, 3)))),omega(1001:2000), -real(fftshift(fft(observables(3001:4000, 2)+1i*observables(3001:4000, 3)))))
