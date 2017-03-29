function out = zfun( xel, c, u, l, thetal, sys )
% z function according to (12) in SAGE paper
%
% xel: estimated received corrupted signal along 'l'th path
% c: c function specified by (1)
% u: trasmtted signal
% l: specific path index
% thetal: (estimated) parameters of the 'l'th path

fs = sys.fs;
I = sys.I;
Ta = sys.Ta;    Tf = sys.Tf;

lena = round(Ta*fs);
lenf = round(Tf*fs);

delay = round(thetal.tau*fs);
ud = zeros(length(u),1);
ud(delay+1:end) = u(1:end-delay);

t = (1:1/fs:length(u)/fs)';
for i = 1:I
    temp = conj(ud).*exp(-1i*2*pi*thetal.fdopp*t).*c.'*xel;
    tempsum = sum(temp((i-1)*lenf+1:(i-1)*lenf+lena));
    tempsumm = tempsumm + tempsum;
end

out = tempsumm;

end

