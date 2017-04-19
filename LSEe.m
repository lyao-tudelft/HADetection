function y = LSEe( in, pilot )
% Least Sqaure Error Equalizer

eL = 3;    % Length of the qualizer
delay = 3;  % Delay <= eL
p = length(in) - delay;

R = toeplitz(in(eL+1:p), in(eL+1:-1:1));
S = pilot(eL+1-delay:p-delay);
size(R);

f = inv(R'*R)*R'*S;
y = filter(f,1,in);
% dec = sign(y);

% figure;
% scatter(real(dec),imag(dec), 'ro'); hold on;
% scatter(real(pilot), imag(pilot),'+');

% error = sum(abs(dec(delay+1:end)-pilot(1:end-delay)));

end

