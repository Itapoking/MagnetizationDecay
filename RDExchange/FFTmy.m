function [freq, spec] = FFTMy(t, signal, phase0, phase1)
    dt = t(2);
    sig =signal;
    N  = size(sig,1);  
    vNyq = 1/(2*dt);      
    zfFactor = 8;
    N_zf = zfFactor * N;        
    sig = [sig; zeros(N_zf-N,1)];
    dv = 1/(N_zf*dt);          
    freq = -vNyq + (0:N_zf-1).' * dv;
    spec = fftshift(fft(sig));
    size(spec)
    size(exp(-1i * (phase0 + phase1 * (0: N_zf-1)') * pi / 180))
    spec = spec.*exp(-1i * (phase0 + phase1 * (0: N_zf-1)') * pi / 180); 
end

