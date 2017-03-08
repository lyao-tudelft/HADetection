function plotTime( in, fs )
% This function plots the input signal "in" in timescale.
% Used for debugging.

len = length(in);
dt = 1/fs/(1e-6);
tmax = (len-1)/fs/(1e-6);
t = linspace(0, tmax, len);

figure();
plot(t, in);
xlabel('Time in us'); ylabel('Magnitued');

end

