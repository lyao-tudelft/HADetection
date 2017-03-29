function out = s( theta, u, sys )
% Received signal without noise corruption along path specified by theta
%
% theta: parameters of specified path
% u: transmitted signal

M = sys.M;
d = sys.d;
fs = sys.fs;

ud = zeros(length(u),1);
delay = round(theta.tau*fs);
ud(delay+1:end) = u(1:end-delay);   % Delayed version of transmitted sinal due to path delay

t = (1:1/fs:length(u)/fs)';
out = zeros(M, length(u));
c = zeros(M,1);
for i = 1:M
    c(i) = exp(1i*2*pi*(i-1)*d*cos(theta.phi));
    out(i,:) = c(i)*theta.amp*exp(1i*2*pi*theta.fdopp.*t).*ud;
end

end

