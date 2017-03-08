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
delayVector = (0.05:0.05:0.2)*1e-6;  % Discrete delays of four-path channel (s)
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

%% Setup Rician channel model
ricChan1 = comm.RicianChannel( ...
    'SampleRate',              sampleRate50MHz, ...
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
    'SampleRate',              sampleRate50MHz, ...
    'PathDelays',              delayVector_a, ...
    'AveragePathGains',        gainVector_a, ...
    'KFactor',                 KFactor_a, ...
    'DirectPathDopplerShift',  specDopplerShift_a, ...
    'MaximumDopplerShift',     maxDopplerShift_a, ...
    'RandomStream',            'mt19937ar with seed', ...
    'Seed',                    100, ...
    'PathGainsOutputPort',     true, ...
    'Visualization',            'off');

Rs = ricChan1.SampleRate;
%% Create modulator & demodulator
modulator = comm.QPSKModulator( ...
    'BitInput',     true, ...
    'PhaseOffset',  pi/4);

demodulator = comm.QPSKDemodulator( ...
    'PhaseOffset',  pi/4, ...
    'BitOutput',    true);

%% Generate a frame of message with 100bits per frame
% payloadLen = 100; frameLength = payloadLen*2;
frameLength = 100;
rb = ricChan1.SampleRate/2;           % Bit rate
Tb = 1/rb;                             % Duration of per bit
Lb = round(Tb*ricChan1.SampleRate);    % Length of per bit after sampling
msg = zeros(frameLength*Lb, 1);
for k = 1:frameLength
    msg( (k-1)*Lb+1:k*Lb ) = randi([0,1],1);
end

% plotTime(msg, ricChan1.SampleRate);
%% Modulate the signal with QPSK scheme with carrier frequency of 240MHz
signal_bb = modulator(msg);        % Baseband modulated signal

fc = 24e6;
t = 0:1/ricChan1.SampleRate:(length(signal_bb)-1)/ricChan1.SampleRate;
signal  = real(signal_bb).*cos(2*pi*fc*t') - 1i.*imag(signal_bb).*sin(2*pi*fc*t');   % Signal has been transated to carrier frequency

% figure('Name','Message and Signal');
% subplot(3,1,1);
% plotSpectrum(msg,Rs);
% subplot(3,1,2);
% plotSpectrum(signal, Rs);
% subplot(3,1,3);
% plotSpectrum(signal_bb, Rs);

%% Pass the signal into channels
[signal_out1,~] = ricChan1(signal);     % Signal from Tx antenna 1
sdr_out1 = sdr(signal_out1, fc, ricChan1.SampleRate);
msg_out1 = demodulator(sdr_out1);

% figure('Name','OUT Message, Signal, and SDR out signal');
% subplot(3,1,2);
% plotSpectrum(signal_out1,Rs);
% subplot(3,1,1);
% plotSpectrum(msg_out1, Rs);
% subplot(3,1,3);
% plotSpectrum(sdr_out1, Rs);

[signal_out2, ~] = ricChan2(signal);    % Signal from Tx antenna 2
sdr_out2 = sdr(signal_out1, fc, ricChan1.SampleRate);
msg_out2 = demodulator(sdr_out2);

%% Setup a timescope system object to view signal magnitude
signalScope = dsp.TimeScope( ...
    'SampleRate', ricChan1.SampleRate, ...
    'TimeSpan',   frameLength/2/ricChan1.SampleRate, ... % One frame span
    'Name',       'Multipath Gain', ...
    'ShowGrid',   true, ...
    'YLimits',    [-5 5], ...
    'YLabel',     'Signal Magnitude', ...
    'MaximizeAxes','On');

% signalScope([msg_out1, msg]);

% Intiail Nulling begins
% h1 = signal_out1./signal;           % Channel estimation for channel 1
% h2 = signal_out2./signal;           % Channel estimation for channel 2

%% Channel estimation based on Least Square Estimation
% Channel estimation
% Training vector
vectorTrain = zeros(2, length(signal));
for k = 1:length(signal)
    vectorTrain(1,k) = signal(k);
    vectorTrain(2,k) = signal(k);
end

% Received trainging vector
vectorReceive = zeros(1, length(signal));
for p = 1:length(signal)
    vectorReceive(1,p) = signal_out1(p) + signal_out2(p);
end

% Least square estimation for CSI
% H = vectorReceive*vectorTrain'*inv(vectorTrain*vectorTrain');

%% Pre-coding
% % msg_n = -h1./h2.*msg;         % Pre-coding
% signal_n = -h1./h2.*signal;
% 
% signal_out1 = ricChan1(signal);     % Antenna 1 transmits the same signal
% signal_out2 = ricChan2(signal_n);   % Antenna 2 ransmits pre-coded signal
% msg_n = demodulator(signal_out1+signal_out2);  % 
% 
% figure(); plot(signal_out1+signal_out2);
% signalScope(signal_out1+signal_out2);