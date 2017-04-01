%% Project simulation.
% This file simulates all the procedure that a message frame may go through
% in the project, including its generation, baseband modulation, passing 
% through a multipath fading channel, reception in SDR dongle,
% demodulation, channel estimation and nullifying, while channel estimation
% and nullifying remains unfinished.

%% Configurations for Rician multi-path channel
sampleRate50MHz = 50e6;              % Sample rate of 50M Hz
sampleRate500KHz = 500e3;            % Sample rate of 500K Hz
sampleRate20KHz  = 20e3;             % Sample rate of 20K Hz
maxDopplerShift  = 200;              % Maximum Doppler shift of diffuse components (Hz)
% delayVector = (0.05:0.05:0.2)*1e-6;  % Discrete delays of four-path channel (s)
gainVector  = [-0.5 -1 -1.5 -2];     % Average path gains (dB)
% One millisecond of delay corresponds to about 300 meters of distance
% difference. For an indoor enironment, path length is no more than dozens
% of meters, which means less than 1 millisecond of delay.

KFactor = 10;            % Linear ratio of specular power to diffuse power
specDopplerShift = 100;  % Doppler shift of specular component (Hz)

% Below define the channel for another Tx antenna
maxDopplerShift_a  = 200;      % Maximum Doppler shift of diffuse components (Hz)
delayVector_a = (0:0.04:0.12)*1e-6; % Discrete delays of four-path channel (s)
gainVector_a  = [0 -0.45 -0.9 -1.35];  % Average path gains (dB)
KFactor_a = 9;            % Linear ratio of specular power to diffuse power
specDopplerShift_a = 98;  % Doppler shift of specular component (Hz)

fs = 3.2e6;
% delayVector = [1/fs,3/fs,5/fs,7/fs];
delayVector = [2/fs 5/fs];
gainVector  = [-2 -3];
%% Setup Rician channel model
ricChan1 = comm.RicianChannel( ...
    'SampleRate',              fs, ...
    'PathDelays',              delayVector, ...
    'AveragePathGains',        gainVector, ...
    'KFactor',                 KFactor, ...
    'DirectPathDopplerShift',  specDopplerShift, ...
    'MaximumDopplerShift',     maxDopplerShift, ...
    'RandomStream',            'mt19937ar with seed', ...
    'Seed',                    100, ...
    'PathGainsOutputPort',     true, ...
    'Visualization',            'off');

ricChan2 = comm.RicianChannel( ...
    'SampleRate',              fs, ...
    'PathDelays',              delayVector_a, ...
    'AveragePathGains',        gainVector_a, ...
    'KFactor',                 KFactor_a, ...
    'DirectPathDopplerShift',  specDopplerShift_a, ...
    'MaximumDopplerShift',     maxDopplerShift_a, ...
    'RandomStream',            'mt19937ar with seed', ...
    'Seed',                    100, ...
    'PathGainsOutputPort',     true, ...
    'Visualization',            'off');

%% Create modulator & demodulator
modulator = comm.QPSKModulator( ...
    'BitInput',     true, ...
    'PhaseOffset',  pi/4);

demodulator = comm.QPSKDemodulator( ...
    'PhaseOffset',  pi/4, ...
    'BitOutput',    true);

%% Generate a sequency of message who contains nFrame frames 
% payloadLen = 100; frameLength = payloadLen*2;
frameLength = 100;
nFrame = 1000;
msgLength = nFrame*frameLength;
rb = fs;           % Bit rate
Tb = 1/rb;                             % Duration of per bit
msg_dur = msgLength*Tb;
Lb = round(Tb*fs);    % Length of per bit after sampling
frame = zeros(frameLength*Lb, 1);
msg = zeros(msgLength*Lb, 1);
for k = 1:frameLength
    frame( (k-1)*Lb+1:k*Lb ) = randi([0,1],1);
end
msg = repmat(frame, nFrame, 1);

% constellation(modulator);
% plotTime(msg, fs);
%% Modulate the signal with QPSK scheme with carrier frequency of 240MHz
signal_bb = modulator(msg);        % Baseband modulated signal

fc = 24e6;
t = 0:1/fs:(length(signal_bb)-1)/fs;
% signal  = real(signal_bb).*cos(2*pi*fc*t') - 1i.*imag(signal_bb).*sin(2*pi*fc*t');   % Signal has been transated to carrier frequency
signal = signal_bb.*cos(2*pi*fc*t');

% figure('Name','Message and Signal');
% subplot(3,1,1);
% plotSpectrum(msg,fs);
% subplot(3,1,2);
% plotSpectrum(signal, fs);
% subplot(3,1,3);
% plotSpectrum(signal_bb, fs);

%% Pass the signal into channels
[signal_out1,~] = ricChan1(signal);     % Signal from Tx antenna 1
sdr_out1 = sdr(signal_out1, fc, fs);
msg_out1 = demodulator(sdr_out1);

[msg_bbout,~] = ricChan1(msg);
signal_outn1 = msg_bbout + randn(length(msg_bbout),1);

% figure('Name','OUT Message, Signal, and SDR out signal');
% subplot(3,1,2);
% plotSpectrum(signal_out1,fs);
% subplot(3,1,1);
% plotSpectrum(msg_out1, fs);
% subplot(3,1,3);
% plotSpectrum(sdr_out1, fs);

[signal_out2, ~] = ricChan2(signal);    % Signal from Tx antenna 2
sdr_out2 = sdr(signal_out2, fc, fs);
msg_out2 = demodulator(sdr_out2);

%% Equalize
dfeObj = dfe(5,3,lms(0.01));
% Set the signal constellation
dfeObj.SigConst = modulator(msg)';
% Maintain continuity between calls to equalize
dfeObj.ResetBeforeFiltering = 0;
% Define initial coefficients to help convergence
dfeObj.Weights = [0 1 0 0 0 0 0 0];

eqRxSig = equalize(dfeObj,signal_out1);

% plot(real(eqRxSig),imag(eqRxSig),'+');


%% Setup a timescope system object to view signal magnitude
signalScope = dsp.TimeScope( ...
    'SampleRate', fs, ...
    'TimeSpan',   frameLength/2/fs, ... % One frame span
    'Name',       'Multipath Gain', ...
    'ShowGrid',   true, ...
    'YLimits',    [-5 5], ...
    'YLabel',     'Signal Magnitude', ...
    'MaximizeAxes','On');

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