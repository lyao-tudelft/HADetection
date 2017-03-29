function out = xel( y, l, theta, u, sys )
% Estimated received signal that's corrupted by White Gaussian noise along
% the 'l'th path . See SAGE paper (13)

% y: received signal with noise corruption after superposition
% theta: parameters of all paths
% u: transmitted signal

fs = sys.fs;
M = sys.M;
d = sys.d;

beta = 1;
L = length(theta);

% out = zeros(M, length(u));
sig = cell(L, 1);
summation = zeros(M, length(u));
for i = 1:L
    sig{i} = s( theta(i), u, sys );
    summation = summation + sig{i};
end

out = sig{l} + beta*(y - summation);


end

