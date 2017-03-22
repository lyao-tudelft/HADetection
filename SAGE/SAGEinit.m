function [ output_args ] = SAGEinit( y, L, M, I, u, fs, d, Ta, Tf )
% Initialization phase of SAGE algorithm
%
% y: superimposed received noise corrupted signal
% L: estimated number of paths
% M: number of receiving antennas
% I: number of observing windows
% u: transmitted signal
% fs: sample rate
% M: number of receving antennas
% d: multiple to wavelength denating the distance among receiving antennas
% Ta: observatin interval duration
% Tf: interval of succesive observation intervals

lena = round(Ta*fs);
lenf = round(Tf*fs);

zero = cell(L,1);
for i = 1:L
    zero{i} = 0;
end

theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', zero);   % Initialize estimated parameters


for miu = -(L-1):0
    l = M+miu;

    xel_last = xel( y, l, theta(l), u, fs, M, d );     % Estimated x with estimated paramters from the last iteration

    %% Estimate tau
    taumax = 2e-6;
    tau = 0:20:taumax;
    %temps = zeros(length(u),1);
    tempsum = zeros(length(tau),1);
    tempsumm = zeros(length(tau),1);
    cou = 1;
    ud = zeros(length(u),1);
    
    % Compute terms within the braces in (16)
    for delay = round(tau*fs)
        ud(delay+1:end) = u(1:end-delay);   % Delayed version of transmitted sinal due to path delay
        
        for m = 1:M
            temps = conj(ud).*xel_last(m,:);
            for i = 1:I
                %temp = (sum(conj(ud((i-1)*lenf+1:(i-1)*lenf+1+lena)).*xel_last(m,(i-1)*lenf+1:(i-1)*lenf+1+lena)))^2;
                temp = (sum(temps((i-1)*lenf+1:(i-1)*lenf+1+lena)))^2;
                tempsum(cou) = tempsum(cou)+temp;
            end
            tempsumm(cou) = tempsumm(cou) + tempsum(cou);
        end
        cou = cou + 1;
    end
    
    % Find the tau that maximize it
    [~,Itau] = max(tempsumm);
    tau_e = 0+(Itau-1)*(taumax/20);    % The delay estimation
    theta(l).tau = tau_e;
    
    %% Estimate phi
    delay_a = round(theta(l).tau*fs);
    ud(delay_a+1:end) = u(1:end-delay_a);   % Delayed version of transmitted signal with ...
                                            % estimated delay from the last iteration
    
    phimax = 2*pi;
    phi = 0:90:phimax;
    c = zeros(M,length(phi));
    xel_last = xel( y, l, theta(l), u, fs, M, d );
    sumphim = zeros(length(phi),1);
    
    % Computing the terms within the braces in (17) for each phi
    for j = 1:length(phi)
        for i = 1:M
            c(i,j) = exp(1i*2*pi*(i-1)*d*cos(phi(j)));
        end
        
        sumphi = conj(ud).*c(:,j).'*xel_last;
        for i = 1:I
            tempphi = (sum(sumphi((i-1)*lenf+1:(i-1)*lenf+lena)))^2;
            sumphim(j) = sumphim(j)+tempphi;
        end
    end
    
    % Find the phi that maximize it
    [~, Iphi] = max(sumphim);
    phi_e = 0+(Iphi-1)*(phimax/90);     % The phi estimation
    theta(l).phi = phi_e;
    
    
end

end

