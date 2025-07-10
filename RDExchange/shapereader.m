function [t, Bampx, Bampy, points, amplitude, phase] = shapereader(filename, tpulse, power, puloffset)
    datapulse = xlsread(filename);
    points = datapulse(:,1); amplitude = datapulse(:,2); phase = datapulse(:,3);
    t = tpulse/points(end)*points;
    Bampx = power*amplitude.*cos(phase/360*2*pi)/max(amplitude);
    Bampy = power*amplitude.*sin(phase/360*2*pi)/max(amplitude);
end

