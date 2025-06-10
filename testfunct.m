 gridsize = [32, 32, 64];
 boxsize = [16, 16, 16];
 Minit = [0, 0, 1];
 shapesize = [5, 0.5];
 DFTset = [200, 200, 200];


 [Minit, K, mask,x,y,z] = Mtest(gridsize, boxsize, Minit, shapesize, DFTset);
 [Bdipx, Bdipy, Bdipz] = Bdip(Minit, K, mask, DFTset, gridsize);

 Nslice=gridsize(3)/2
 Nstepshow=2
 Nstepz=2
surf(x(:,:,Nslice), y(:,:,Nslice),  -Bdipz(:,:,Nslice));
% quiver3(x(1:Nstepshow:end),y(1:Nstepshow:end),z(1:Nstepz:end),Bdipx(1:Nstepshow:end),Bdipy(1:Nstepshow:end),Bdipz(1:Nstepz:end))