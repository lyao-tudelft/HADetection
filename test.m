f = 10e6;
T = 1/f;
t = linspace(1,5*T,100);
y = cos(2*pi*f*t);

figure();
plotSpectrum(y,20*f);

yt = y.*exp(1i*2*pi*2*f*t);
figure();
plotSpectrum(yt,20*f);