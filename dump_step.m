function dump_step(M, dw, wdip, wrd)
    dump_hdd('Magnetization.bin', M)
    dump_hdd('Offset.bin', dw)
    dump_hdd('Dipolarfrequency.bin', wdip)
    dump_hdd('RDfrequency.bin', wrd)

end