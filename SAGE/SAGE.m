function theta_out = SAGE( y, u, theta, sys )
% Main phase of SAGE algorithm
%
% y: superimposed received signal
% u: transmitted signal
% theta: intialized channel parameters
% sys: system information
%
% theta_out: result of estimated channel parameters

L = sys.L;
d = sys.d;
M = sys.M;
I = sys.I;
Ta = sys.Ta;
fs = sys.fs;
lena = round(Ta*fs);

distance = inf;
threshold = 1e-3;
miu = 0;
miumax = 10;

theta_temp = theta;
% while distance >= threshold && miu <= miumax
while miu < miumax
    
    miu = miu+1
    theta_last = theta_temp;
    for l = 1:L
        xel_last = xel( y, l, theta, u, sys );
        
        % Estimate time delay tau
        c = zeros(M,1);
        for m = 1:M
            c(m) = exp(1i*2*pi*(m-1)*d*cos(theta_temp(l).phi));
        end
        
        taumax = 4e-6;
        tau = linspace(0,taumax,500);
        theta_t = theta_temp;
        for i = 1:length(tau)
            theta_t(l).tau = tau(i);
            z(i) = abs(zfun( xel_last, c, u, theta_t(l), sys ));
        end
        [~,Itau] = max(z);
        tau_e = tau(Itau);
        theta_temp(l).tau = tau_e;
        
        % Estimate DoA phi
        phimax = 2*pi;
        phi = linspace(0,phimax,5);
        
        cphi = zeros(M,length(phi));
        theta_t = theta_temp;
        for j = 1:length(phi)
            theta_t(l).phi = phi(j);
            for m = 1:M
                cphi(m,j) = exp(1i*2*pi*(m-1)*d*cos(theta_t(l).phi));
            end
            
            z(j) = abs(zfun( xel_last, cphi(:,j), u, theta_t(l), sys ));
        end
        [~, Iphi] = max(z);
        phi_e = phi(Iphi);
        theta_temp(l).phi = phi_e;
        
        % Estimate doppler frequency fdopp
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
        [~, Ifd] = max(z);
        fd_e = fd(Ifd);
        theta_temp(l).fdopp = fd_e;
        
        % Estimate amplitude amp
        Pu = sum(u.^2)/length(u);
        a_e = 1/(I*norm(cfd)^2*lena*Pu)*zfun(xel_last, cfd, u, theta_temp(l), sys);
        theta_temp(l).amp = a_e;
        
    end
    theta_out = theta_temp;
    
    for i = 1:L
        distance = sum(abs(theta_out(l).tau - theta_last(l).tau)) + sum(abs(theta_out(l).amp - theta_last(l).amp));
    end
end
    
end

