function Mfull = offsetoptim(M0, M00,  T1, T2, T, dw, time, npoints, k)
    ndw = size(dw, 1);
    Mfull = zeros(npoints, 3*2*ndw+1);


    
    for i=1:ndw
        dwoff = dw(i,:)
        [t, Mtraj] = ode113(@(t, M) twodensBMConnell(M, dwoff, T1, T2, T, 0, 0, M00, k), linspace(0, time, npoints), M0);
        Mfull(:, 2 + 6*(i-1)) = Mtraj(:,1);
        Mfull(:, 3 + 6*(i-1)) = Mtraj(:,2);
        Mfull(:, 4 + 6*(i-1)) = Mtraj(:,3);
        Mfull(:, 5 + 6*(i-1)) = Mtraj(:,4);
        Mfull(:, 6 + 6*(i-1)) = Mtraj(:,5);
        Mfull(:, 7 + 6*(i-1)) = Mtraj(:,6);

    end

    Mfull(:, 1) = t;
    
end
