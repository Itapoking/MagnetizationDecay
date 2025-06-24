function dump_hdd(filename, array)
    fid = fopen(filename, 'a');
    fwrite(fid, array, 'double');
    fclose(fid);
end