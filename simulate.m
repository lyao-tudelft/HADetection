M = 4; % Alphabet size for modulation
msg = randi([0 M-1],2500,1); % Random message
hMod = comm.QPSKModulator('PhaseOffset',0);
modmsg = hMod(msg); % Modulate using QPSK
chan = [.986; .845; .237; .123+.31i]; % Channel coefficients
filtmsg = filter(chan,1,modmsg); % Introduce channel distortion

dfeObj = dfe(5,3,lms(0.01));
% Set the signal constellation
dfeObj.SigConst = hMod((0:M-1)')';
% Maintain continuity between calls to equalize
dfeObj.ResetBeforeFiltering = 0;
% Define initial coefficients to help convergence
dfeObj.Weights = [0 1 0 0 0 0 0 0];

eqRxSig = equalize(dfeObj,filtmsg);

initial = eqRxSig(1:200);
plot(real(initial),imag(initial),'+')
hold on
final = eqRxSig(end-200:end);
plot(real(final),imag(final),'ro')
legend('initial', 'final')

y = equalizer(10,filtmsg,modmsg);

% figure;
% initialy = y(1:200);
% plot(real(initialy),imag(initialy),'+')
% hold on
% finaly = y(end-200:end);
% plot(real(finaly),imag(finaly),'ro')
% legend('initial', 'final')