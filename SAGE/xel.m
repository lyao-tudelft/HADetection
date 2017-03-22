function out = xel( y, l, theta, u, fs, M, d )
% Estimated received signal that's corrupted by White Gaussian noise along
% the 'l'th path . See SAGE paper (13)

% y: received signal with noise corruption after superposition
% theta: parameters of all paths
% u: transmitted signal
% fs: sample rate
% M: number of receving antennas
% d: multiple to wavelength denating the distance among receiving antennas

beta = 1;
L = length(theta);

% out = zeros(M, length(u));
sig = cell(L, 1);
summation = zeros(M, length(u));
for i = 1:L
    sig{i} = s( theta(i), u, fs, M, d );
    summation = summation + sig{i};
end

out = sig{l} + beta*(y - summation);


end

