function dMdt = twodensBMConnell(M, dw, T1, T2, Trd, B1x, B1y, M00, k)
    Max = M(1); May = M(2); Maz = M(3);
    Mbx = M(4); Mby = M(5); Mbz = M(6);
    dwa = dw(1); dwb = dw(2); M00a = M00(1); M00b = M00(2);
    T1a = T1(1); T1b = T1(2); T2a = T2(1); T2b = T2(2);
    lrd =  1/Trd;

    dMdt = [ -Max/T2a + dwa*May + B1y*Maz - lrd*Maz*(Max + Mbx) - k*Max + k*Mbx;
             -May/T2a - dwa*Max - B1x*Maz  - lrd*Maz*(May + Mby) - k*May + k*Mby;
             -B1y*Max + B1x*May - (Maz - M00a)/T1a + lrd*(Max^2 + May^2 + Max*Mbx + May*Mby) - k*Maz + k*Mbz;
             -Mbx/T2b + dwb*Mby + B1y*Mbz - lrd*Mbz*(Max + Mbx) - k*Mbx + k*Max;
             -Mby/T2b - dwb*Mbx - B1x*Mbz  - lrd*Mbz*(May + Mby) - k*Mby + k*May;
             -B1y*Mbx + B1x*Mby - (Mbz - M00b)/T1b + lrd*(Mbx^2 + Mby^2 + Max*Mbx + May*Mby) - k*Mbz + k*Maz;];

end