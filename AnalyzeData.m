function [Mx, My, Mz, Mabs, wrdx, wrdy, wrdz, t] = AnalyzeData(M, wrd, t)
    [Mx, My, Mz, Mabs] = M_averaging(M);
    [wrdx, wrdy, wrdz] = M_averaging(wrd);
    dump_average(M, dw, wdip, wrd)
    % [wdipx, wdipy, wdipz, wdipabs] = M_averaging(wdip);
end

