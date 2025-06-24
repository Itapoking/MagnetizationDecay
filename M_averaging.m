function [Mx, My, Mz, Mnorm] = M_averaging(M, mask)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  M = [Mx, My, Mz] - 3 3D arrays with magnetization values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = reshape(M, [size(mask) 3]);
    sizemask = size(mask);
    Vgrid = sizemask(1)*sizemask(2)*sizemask(3);
    Mx = mean(M(:,:,:,1), 'all')*Vgrid/sum(mask(:,:,:),"All"); 
    My = mean(M(:,:,:,2), 'all')*Vgrid/sum(mask(:,:,:),"All"); 
    Mz = mean(M(:,:,:,3), 'all')*Vgrid/sum(mask(:,:,:),"All"); 
    Mnorm = sqrt(Mx^2 + My^2 + Mz^2);
    
end

