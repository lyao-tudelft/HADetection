fs = 3.2e6;

numAnt = 'M';       M = 2;      % Number of receiving antennas
numWin = 'I';       I = 20;     % Number of observing windows
sampFreq = 'fs';    rs = fs;   % Sample Rate
lenWin = 'Ta';      Ta = 0.001;    % Length of observing windows
intervWin = 'Tf';   Tf = 0.0015;    % Interval of observing windows
dAnt = 'd';         d = 1;      % Multiple of wavelength denoting the distance among antennas
numPath = 'L';      L = 4;

sys = struct(numAnt, M, numWin, I, sampFreq, rs, lenWin, Ta, intervWin, Tf, dAnt, d, numPath, L);

u = msg;

delay = [2/fs 5/fs];
amp = [0.2 0.5];
DoA = [pi/3 pi/6];
fd = [0 0];

CSI = struct('delay',delay,'amp',amp,'DoA',DoA,'fd',fd);

s = signalGenerator( u, CSI, sys );

zero = cell(L,1);
one = cell(L,1);
for i = 1:L
    zero{i} = 0;
    one{i} = 0;
end
theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', one);
theta = SAGE( s, msg, theta, sys );