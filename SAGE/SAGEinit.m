function [ output_args ] = SAGEinit( y, L, u, sys )
% Initialization phase of SAGE algorithm. See details in SAGE paper IV.B
%
% y: superimposed received noise corrupted signal
% L: estimated number of paths
% u: transmitted signal
% sys: system information

Ta = sys.Ta;    Tf = sys.Tf;
fs = sys.fs;
M = sys.M;
d = sys.d;
I = sys.I;

lena = round(Ta*fs);
lenf = round(Tf*fs);

zero = cell(L,1);
for i = 1:L
    zero{i} = 0;
end

theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', zero);   % Initialize estimated parameters


for miu = -(L-1):0
    l = M+miu;

    xel_last = xel( y, l, theta(l), u, sys );     % Estimated x with estimated paramters from the last iteration

    %% Estimate tau
    taumax = 2e-6;
    tau = linspace(0,taumax,20);
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
    
    %% Estimate phi
    delay_a = round(theta(l).tau*fs);
    ud(delay_a+1:end) = u(1:end-delay_a);   % Delayed version of transmitted signal with ...
                                            % estimated delay from the last iteration
    
    phimax = 2*pi;
    phi = linspace(0,phimax,90);
    c = zeros(M,length(phi));
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
    
    %% Estimate doppler frequency
    fdmax = 10;
    fd = linspace(0,fdmax,10);  % Doppler frequency
    z = zfun( xel, c, u, l, thetal, sys );
end

end

