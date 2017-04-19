%% Simulate multipath channel and conduct estimation with SAGE and Wiener Filter
fs = 3.2e6;

%% Generate PN sequence awaiting to be sent
frameLength = 2000;
nFrame = 1;
msgLength = nFrame*frameLength;

rb = fs;           % Bit rate
Tb = 1/rb;                             % Duration of per bit
Lb = round(Tb*fs);    % Length of per bit after sampling
Tmsg = Tb*msgLength;

H = comm.PNSequence('MaximumOutputSize', [frameLength 1],...
                    'VariableSizeOutput', 1, ...
                    'Polynomial', [12 6 4 1 0], ...
                    'InitialConditions', [1 0 0 0 0 0 1 0 0 0 0 0]);
PN = 2*H(frameLength)-1;
% PN = H(frameLength);

frame = zeros(frameLength*Lb, 1);
msg = zeros(msgLength*Lb, 1);
for k = 1:frameLength
    frame( (k-1)*Lb+1:k*Lb ) = PN(k)*ones(Lb,1);
end
msg = repmat(frame, nFrame, 1);
% msg(1) = 1;

%% System Information
numAnt = 'M';       M = 2;      % Number of receiving antennas
numWin = 'I';       I = 20;     % Number of observing windows
sampFreq = 'fs';    rs = fs;   % Sample Rate
lenWin = 'Ta';      Ta = Tmsg/(I+1);    % Length of observing windows
intervWin = 'Tf';   Tf = Tmsg/I;    % Interval of observing windows
dAnt = 'd';         d = 1;      % Multiple of wavelength denoting the distance among antennas
numPath = 'L';      L = 4;

sys = struct(numAnt, M, numWin, I, sampFreq, rs, lenWin, Ta, intervWin, Tf, dAnt, d, numPath, L);

u = msg;

%% Specify Channel State Information
delay = [50/fs 100/fs 170/fs];
amp = [0.2 0.5 0.3];
DoA = [0 0 0];
fd = [0 0 0];

CSI = struct('delay',delay,'amp',amp,'DoA',DoA,'fd',fd);

%% Disperse and attenuate the signal in the channel
[s, sn] = signalGenerator( u, CSI, sys );

%% SAGE
% zero = cell(L,1);
% one = cell(L,1);
% for i = 1:L
%     zero{i} = 0;
%     one{i} = 0;
% end
% theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', one);

% theta = SAGEinit( s, msg, sys );
% theta = SAGE( s, msg, theta, sys );

%% Wiener Filter

% Auto-correlation of sent signal
[rx, lagrx] = xcorr(msg, 'unbiased');
rx = rx(lagrx>=0);
lagrx = lagrx(lagrx>=0);

% Cross-correlation between sent signal and received signal
[rdx, lagrdx] = xcorr(s(:,1), msg, 'unbiased');
rdx = rdx(lagrdx>=0);
lagrdx = lagrdx(lagrdx>=0);

Px = fft(rx(1:round(length(rx)*0.7)));
Pdx = fft(rdx(1:round(length(rdx)*0.7)));
H = Pdx./Px;
h = ifft(H);    % Estimated channel impulse response in Minimum Least Square Sense

%% Plot result
delayt = zeros(length(h),1);
delayt(delay*fs) = amp;

figure
plot([1:length(h)]/fs,h);
hold on;    grid on;
title('Channel Impulse Response by Wiener');
xlabel('Time / s');
ylabel('Amplitude');
plot([1:length(h)]/fs,delayt);
axis([0 8e-5 -0.1 0.7]);
legend('Estimated','Real');
% axis([0 length(h)*0.2/fs -0.1 1]);
% plot(delay, amp);

% figure
% plot(rx(1:round(length(rx)*0.7)));
% hold on;
% plot(rdx(1:round(length(rdx)*0.7)));