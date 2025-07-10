clear all; clc;

phi = 6/360*2*pi; theta = 0/360*2*pi;
phi2 = 6/360*2*pi; theta2 = 90/360*2*pi;
enh = [500, 10];
M00 = [1, 1];
M0 = [enh(1)*sin(phi)*cos(theta), enh(1)*sin(phi)*sin(theta), enh(1)*cos(phi), enh(2)*sin(phi2)*cos(theta2), enh(2)*sin(phi2)*sin(theta2), enh(2)*cos(phi2)];

k = 800;

dw =  2*pi*500*([0, -0.2;]);

T1 = [15, 15]; T2 = 0.1*[1, 1];
T = 0.085;

% a = './reburp1000.xlsx';
% tsoft = 3237*1e-6;  % 12 ms full inversion
% Bsoft = 1896*2*pi;
% puloffset = 0*500;
% [tpulse, Bampx, Bampy, points, amplitude, phase] = shapereader(a, tsoft, Bsoft, puloffset);
% % 
% % 
npotinssim=1e5;
% Mfull = offsetoptimreburp(M0, M00,  T1, T2, T, dw, 0.4, npotinssim, Bampx,Bampy, tpulse, k);

Mfull = offsetoptim(M0, M00,  T1, T2, T, dw, 0.4, npotinssim, k);

t = Mfull(:,1);
signal = zeros(npotinssim, size(dw,1));
specfull = zeros(8*npotinssim, size(dw,1));

phase0 = 0*[10];
for n = 1:size(dw,1)
    signaltemp = Mfull(:, 2 + 6*(n - 1)) + 1i*Mfull(:, 3 + 6*(n - 1)) +  Mfull(:, 5 + 6*(n - 1)) + 1i*Mfull(:, 6 + 6*(n - 1));
    signal(:,n) = signaltemp;
    [freq, spec] = FFTmy(t, signaltemp, phase0(n), 0.0);
    specfull(:,n) = spec;
end
chemfreq = freq/500;


plot(chemfreq, specfull(:,:), 'LineWidth', 1);
ylabel('Intensity, a.u.','FontSize',20);
xlabel('Chem. shift, ppm','FontSize',20);
xlim([-10 70])
% lgd = legend({'-0.1 ppm', '-2 ppm','-4 ppm','-7 ppm', '-15 ppm','-24 ppm','-60 ppm'}, 'Location', 'NorthWest'); drawnow;
% set(lgd, 'FontSize', 12)
% legend('Location', 'southeast')
% newcolors = {'red', 'black', 'green', 'blue', 'magenta', 'yellow', 'magenta'};
% colororder(newcolors)



% ylim([-2e3 3e3])
% plot(t, real(signal(:,1)), 'LineWidth',2)
% plot(t, sqrt(Mtraj(:,1).^2+Mtraj(:,2).^2+Mtraj(:,3).^2)+sqrt(Mtraj(:,4).^2+Mtraj(:,5).^2+Mtraj(:,6).^2), 'LineWidth',2)
% plot(t, sqrt(Mtraj(:,1).^2+Mtraj(:,4).^2+Mtraj(:,2).^2+Mtraj(:,5).^2), 'LineWidth',2)
% plot(t, sqrt(Mtraj(:,3).^2+Mtraj(:,6).^2), 'LineWidth',2)
