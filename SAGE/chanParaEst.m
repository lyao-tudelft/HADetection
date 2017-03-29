numAnt = 'M';       M = 2;      % Number of receiving antennas
numWin = 'I';       I = 10;     % Number of observing windows
sampFreq = 'fs';    fs = 2e6;   % Sample Rate
lenWin = 'Ta';      Ta = 10;    % Length of observing windows
intervWin = 'Tf';   Tf = 15;    % Interval of observing windows
dAnt = 'd';         d = 1;      % Multiple of wavelength denoting the distance among antennas

sys = struct(numAnt, M, numWin, I, sampFreq, fs, lenWin, Ta, intervWin, Tf, dAnt, d);