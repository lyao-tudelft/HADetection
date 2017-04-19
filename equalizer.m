function  out  = equalizer( L, input, pilot )
% Fractionally Spaced Equalizer according to paper "SAM: Enabling Practical Spatial Multiple Access in
% Wireless LAN" by K. Tan, H. Liu, etc., Figure 8(b).

% Initializing equalizer
c = zeros(2*L+1,1);
c(L+1) = 1;
K = 10;
nSamp = length(input);
b = zeros(2*L+1,1);        % Buffer
bf = b;                  % Flipped buffer
pilote = zeros(nSamp,1);
error = zeros(K, 1);
ss = 0.01;                % Step size

% Training
for i = 1:nSamp
    
    for k = 0:K-1
        
        if i+k-L <= 0
            b(1:-(i+k-L)+1) = 0;
            b(-(i+k-L)+2:end) = input(1:(2*L+1+(i+k-L)-1));
        elseif i+k+L > nSamp
            b(end-(i+k+L-nSamp):end) = 0;
            b(1:end-(i+k+L-nSamp)) = input(i+k-L:end);
        else
            b = input(i+k-L:i+k+L);
        end
        
        if i+k > nSamp
            error(k+1) = 0;
        else
            pilote(i+k) = c.'*flipud(b);
            error(k+1) = pilot(i+k) - pilote(i+k);
        end
        
        delta = zeros(2*L+1,1);
        w = 0;
        
        bf = flipud(b);
        for j = 1:2*L+1
            delta(j) = delta(j) + ss*error(k+1)*conj(bf(j));
            w = w + abs(bf(j))^2;
        end
    end
    
    for j = 1:2*L+1
        c(j) = c(j) - delta(j)/w;
    end
    
end
    
out = filter(c,1,input);
% e = out-pilot;

cons = comm.ConstellationDiagram();

figure;
% plot(real(input),imag(input),'ro'); hold on;
% plot(real(out),imag(out),'+');
% legend('input','out');

% figure;
% plot(real(pilot),imag(pilot),'+');


end



