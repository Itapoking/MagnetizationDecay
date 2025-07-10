clear all; clc;

phi = 30/360*2*pi; theta = 0/360*2*pi;
phi2 = 30/360*2*pi; theta2 = 90/360*2*pi;
enh = [900, 900];
M00 = [1, 1]
M0 = [enh(1)*sin(phi)*cos(theta), enh(1)*sin(phi)*sin(theta), enh(1)*cos(phi), enh(2)*sin(phi2)*cos(theta2), enh(2)*sin(phi2)*sin(theta2), enh(2)*cos(phi2)];

dw =  2*pi*500*([0.0, 1;]);

T1 = 1*[15, 15]; T2 = 0.1*[1, 1];
T = 0.085;
npotinssim=1e6;
Mfull = offsetoptim(M0, M00,  T1, T2, T, dw, 0.001, npotinssim);

t = Mfull(:,1);
signal = zeros(npotinssim, size(dw,1));
specfull = zeros(8*npotinssim, size(dw,1));

phase0 = 0*[120 120 120 120 120 120];
for n = 1:size(dw,1)
    signaltemp = Mfull(:, 2 + 6*(n - 1)) + 1i*Mfull(:, 3 + 6*(n - 1)) +  Mfull(:, 5 + 6*(n - 1)) + 1i*Mfull(:, 6 + 6*(n - 1));
    signal(:,n) = signaltemp;
    [freq, spec] = FFTmy(t, signaltemp, phase0(n), 0.0);
    specfull(:,n) = spec;
end
chemfreq = freq/500;

Manorm = sqrt(Mfull(:,2).^2 + Mfull(:,3).^2 + Mfull(:,4).^2);
Mbnorm = sqrt(Mfull(:,5).^2 + Mfull(:,6).^2)% + Mfull(:,7).^2);
plot(Mbnorm)

% 
% lH =plot3(Mfull(:,2), Mfull(:,3), Mfull(:,4), 'r-', LineWidth=2)
% aH = ancestor(lH,'axes');
% set(aH,'xLim',[-enh(1) enh(1)])
% set(aH,'yLim',[-enh(1) enh(1)])
% set(aH,'zLim',[-enh(1) enh(1)])
% % hold on
% lH1 = plot3(Mfull(:,5), Mfull(:,6), Mfull(:,7), 'g-', LineWidth=2)
% lgd = legend({'T_r_d = Inf', 'T_r_d = 0.0001'},'Location','NorthWest'); drawnow;
% 
% 
% h = animatedline('MaximumNumPoints', 10);
% 

% view(3);
% 
% hold on
% [X,Y,Z]=sphere(enh(1)/enh(1));
% r = 1;
% X2 = X * r;
% Y2 = Y * r;
% Z2 = Z * r;
% 
% 
% surf(X2,Y2,Z2,'EdgeColor','none', 'FaceColor',[0 0 1], 'FaceAlpha',0.1)
% lgd = legend({'sphere'},'Location','NorthWest'); drawnow;
% kylabel('My','FontWeight', 'Bold','FontSize',20); 
% kxlabel('Mx','FontWeight', 'Bold','FontSize',20);
% kzlabel('Mz','FontWeight', 'Bold','FontSize',20);kgrid;
% daspect([1 1 1])
% axis normal
% 
% dstep=1e2
% Mshape = Mfull(1:dstep:end,2:4)/sqrt(enh(1)^2)
% Mshape1 = Mfull(1:dstep:end,5:7)/sqrt(enh(2)^2)
% grid on;
% h = animatedline('MaximumNumPoints',size(Mshape,1), "color", "red",'LineWidth',2);
% hrd = animatedline('MaximumNumPoints', size(Mshape,1), "Color","green",'LineWidth',2);
% hvec = quiver3(0, 0, 0, 0, 0, 0, 'LineWidth', 2, 'Color', 'blue', 'MaxHeadSize', 0.5);
% 
% xvec = quiver3(0, 0, 0, 1, 0, 0, 'LineWidth', 1, 'Color', 'black', 'MaxHeadSize', 0.1);
% yvec = quiver3(0, 0, 0, 0, 1, 0, 'LineWidth', 1, 'Color', 'black', 'MaxHeadSize', 0.1);
% zvec = quiver3(0, 0, 0, 0, 0, 1, 'LineWidth', 1, 'Color', 'black', 'MaxHeadSize', 0.1);
% lgd = legend({'sphere', 'T_r_d = Inf', },'Location','NorthWest'); drawnow;
% 
% view(3)
% % campos([2, 2, 1.5]);     % положение камеры: между X и Y, с высотой
% % camtarget([0, 0, 0]);    % центр, куда смотрим
% % camup([0, 0, 1]);
% % pause(10)
% axis([-1, 1, -1, 1, -1, 1])
% daspect([1 1 1])
% axis normal
% tpause=0.05
% for k = 1:size(Mshape,1)
%     addpoints(h,Mshape(k,1),Mshape(k,2),Mshape(k,3));
%     addpoints(hrd,Mshape1(k,1),Mshape1(k,2),Mshape1(k,3));
%     set(hvec, 'XData', 0, 'YData', 0, 'ZData', 0, ...
%               'UData', Mshape(k,2)+Mshape1(k,2), 'VData', -Mshape(k,1)-Mshape1(k,1), 'WData', 0);
%     pause(tpause)
%     drawnow
% end
