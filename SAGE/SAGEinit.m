function theta_out = SAGEinit( y, u, sys )
% Initialization phase of SAGE algorithm by successive interference
% cancellation. See details in SAGE paper IV.B
%
% y: superimposed received noise corrupted signal
% u: transmitted signal
% sys: system information
%
% theta_out: initialized channel parameters

Ta = sys.Ta;    Tf = sys.Tf;
fs = sys.fs;
M = sys.M;
d = sys.d;
I = sys.I;
L = sys.L;

lena = round(Ta*fs);
lenf = round(Tf*fs);

zero = cell(L,1);
one = cell(L,1);
for i = 1:L
    zero{i} = 0;
    one{i} = 0;
end

theta = struct('tau', zero, 'phi', zero, 'fdopp', zero, 'amp', one);   % Initialize estimated parameters
theta_temp = theta;

for miu = -(L-1):0
    l = L+miu;
    
    xel_last = xel( y, l, theta, u, sys );     % Estimated x with estimated paramters from the last iteration

    %% Estimate tau
    taumax = 4e-6;
    tau = linspace(0,taumax,500);
    %temps = zeros(length(u),1);
    tempsum = zeros(length(tau),1);
    tempsumm = zeros(length(tau),1);
    cou = 1;
    ud = zeros(length(u),1);
    
    % Compute terms within the braces in (16)
    for delay = round(tau*fs)
        ud(delay+1:end) = u(1:end-delay);   % Delayed version of transmitted sinal due to path delay
        
        for m = 1:M

            temps = conj(ud)'.*xel_last(m,:);
            for i = 1:I
                %temp = (sum(conj(ud((i-1)*lenf+1:(i-1)*lenf+1+lena)).*xel_last(m,(i-1)*lenf+1:(i-1)*lenf+1+lena)))^2;
                temp = abs(sum(temps((i-1)*lenf+1:(i-1)*lenf+1+lena)))^2;
                tempsum(cou) = tempsum(cou)+temp;
            end
            tempsumm(cou) = tempsumm(cou) + tempsum(cou);
        end
        cou = cou + 1;
    end
    tempsum;
    % Find the tau that maximize it
    [~,Itau] = max(tempsumm);
    Itau;
    tau_e = tau(Itau);    % The delay estimation
    theta_temp(l).tau = tau_e;
    
    %% Estimate phi
    delay_a = round(tau_e*fs);
    ud(delay_a+1:end) = u(1:end-delay_a);   % Delayed version of transmitted signal with ...
                                            % estimated delay from the last iteration
    
    phimax = 2*pi;
    phi = linspace(0,phimax,5);
    c = zeros(M,length(phi));
    sumphim = zeros(length(phi),1);
    
    % Computing the terms within the braces in (17) for each phi
    for j = 1:length(phi)
        for i = 1:M
            c(i,j) = exp(1i*2*pi*(i-1)*d*cos(phi(j)));
        end

        sumphi = conj(ud)'.*((c(:,j).')*xel_last);
        for i = 1:I
            tempphi = abs(sum(sumphi((i-1)*lenf+1:(i-1)*lenf+lena)))^2;
            sumphim(j) = sumphim(j)+tempphi;
        end
    end
    sumphim;
    % Find the phi that maximize it
    [~, Iphi] = max(sumphim);
    phi_e = phi(Iphi);     % The phi estimation
    theta_temp(l).phi = phi_e;
    
    %% Estimate doppler frequency
    fdmax = 10;
    fd = linspace(0,fdmax,10);  % Doppler frequency
    
    cfd = zeros(M,1);
    for m = 1:M
        cfd(m) = exp(1i*2*pi*(m-1)*d*cos(theta_temp(l).phi));
    end
    
    theta_t = theta_temp;
    z = zeros(length(fd),1);
    for i = 1:length(fd)
        theta_t(1).fdopp = fd(i);
        z(i) = abs(zfun(xel_last, cfd, u, theta_t(l), sys));
    end
    
    z;
    [~, Ifd] = max(z);
    Ifd;
    fd_e = fd(Ifd);
    theta_temp(l).fdopp = fd_e;
    
    %% Estimate Amplitude
    Pu = sum(u.^2)/length(u);
    a_e = 1/(I*norm(c)^2*lena*Pu)*zfun(xel_last, cfd, u, theta_temp(l),sys);
    theta_temp(l).amp = a_e;
    
end

theta_out = theta_temp;

end

