function out = sdr( in,fc,fs )
% This function simulates the procedure that an input signal goes through 
% after received by the front end of SDR dongle, including frequency
% downward converting and low pass filter.
% 
% At last, a signal "out" including two quadrature signals is outputed.

% Frequency translation on input signal
t = linspace(0,(length(in)-1)/fs,length(in));
t = t';
in_t = real(in).*exp(-1i*2*pi*fc*t) + 1i*imag(in).*exp(-1i*2*pi*fc*t);     

% Low pass filter design
Wp = 0.8/2*fs*2*pi/fs/pi;               % Cut-off frequency normalized by Nyquist frequency
Ws = 0.5*fs*2*pi/fs/pi;                 % Stop band corner frequency
Rp = 0.1;                               % Max passband loss in dB
Rs = 30;                                % Stopband attenuation in dB
[n, Wp] = cheb1ord(Wp, Ws, Rp, Rs);     % Chebyshev filter design
[nu, denu] = cheby1(n, Rp, Wp);         % Filter coefficients

filter_out_r = filter(nu, denu, real(in_t));
filter_out_i = filter(nu, denu, imag(in_t));

out = filter_out_r + 1i*filter_out_i;

% figure('Name','Check filter');
% subplot(2,1,1);
% plotSpectrum(in_t, fs);
% subplot(2,1,2);
% plotSpectrum(in, fs);
end

