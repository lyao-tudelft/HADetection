load('channelSim.mat');
%% SAGE
numAnt = 'M';       M = 2;      % Number of receiving antennas
numWin = 'I';       I = 20;     % Number of observing windows
sampFreq = 'fs';    rs = fs;   % Sample Rate
lenWin = 'Ta';      Ta = 0.001;    % Length of observing windows
intervWin = 'Tf';   Tf = 0.0015;    % Interval of observing windows
dAnt = 'd';         d = 1;      % Multiple of wavelength denoting the distance among antennas
numPath = 'L';      L = 4;

sys = struct(numAnt, M, numWin, I, sampFreq, rs, lenWin, Ta, intervWin, Tf, dAnt, d, numPath, L);

% theta = SAGEinit( signal_out1, signal, sys );
zero = cell(L,1);
one = cell(L,1);
for i = 1:L
    zero{i} = 0;
    one{i} = 0;
end

theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', one);   % Initialize estimated parameters

s1 = signal_outn1;
s2 = signal_outn1*exp(1i*2*pi*d*cos(pi/6));
s = [s1 s2];
theta = SAGE( s, msg, theta, sys );