function [s, s_n] = signalGenerator( u, CSI, sys )
% Compute the noise and multipath attenuated signal as the ouput of
% channel, whose CSI is defined in 'CSI'.
%
% u: the transmitted signal.
% CSI: the predefined channel CSI including amplitude(amp), time delay(delay),
%      Doppler frequency(fd) and drection of arrival(DoA) of each path.
% sys: system information
%
% s: signal received in the receiver without noise
% s_n: noise corrupted s

L = length(CSI.delay);
M = sys.M;
d = sys.d;
fs = sys.fs;


t = linspace(0,round((length(u)-1)/fs), length(u))';
s = zeros(length(u), M);
for m = 1:M
    stemp = zeros(length(u), L);
    for l = 1:L
        delayd = round(CSI.delay(l)*fs);
        ud = zeros(length(u),1);
        ud(delayd+1:end) = u(1:end-delayd);
        
        stemp(:,l) = exp(1i*2*pi*(m-1)*d*cos(CSI.DoA(l)))*CSI.amp(l)*exp(1i*2*pi*CSI.fd(l)*t).*ud;
    end
    s(:,m) = sum(stemp,2);
end

s_n = s + sqrt(0.2)*randn(length(s),M);

end

