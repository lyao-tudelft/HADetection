function plotSpectrum( in,fs )
% Thsi function plots single sided frequency spectrum for input signal "in"
% at sampling frequency "fs".
% Used for debugging.

L = length(in);

in_s = fft(in);
P2 = abs(in_s/L);                   % Double sided spectrum
P1 = P2(1:L/2+1);                   % Single sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:(L/2))/L;
plot(f, P1);
xlabel('Frequency/Hz'); ylabel('Magnitude');

end

