clear

%load('theta.mat');  %results from sage saved in file theta.mat

N=4; %number of CSIs in each sample
N1=3; %number of samples in each time slot/window
N_tot=50;  %total number of channel measurement sample
%N_tot=size(theta,2);
group=N_tot-N1;

n=1:N_tot;
g=1:group;
H_ori=zeros(N,length(n));
H_ori_no=zeros(N,length(n));
H=zeros(N,N1,group);
A=zeros(N,N,group-1);
C=zeros(N,N,group-1);
A_no=zeros(N,N,group-1);
C_no=zeros(N,N,group-1);
lamda_a=zeros(1,group);

%%generate unrelated matrixs to simulate huamn movements
for ii=1:N_tot
    H_ori(:,ii)=rand(N,1)+i*rand(N,1);  %generate 10 random channel measurement sampel
end
%H_ori=theta;  %copy the SAGE result to H

for i=1:group
    m=i;
    for j=1:N1
    H(:,j,i)=H_ori(:,m);  %put N1 samples in each time slot
    m=m+1;
    end
end

%calculate the correlation matrix A and C for amplitudes and phase of the CSI measurements
for ii=1:group
    for i=1:N
        for j=1:N
                temp=corrcoef(real(H(i,:,ii)),real(H(j,:,ii)));
                A(i,j,ii)=temp(1,2);
                temp=corrcoef(imag(H(i,:,ii)),imag(H(j,:,ii)));  
                C(i,j,ii)=temp(1,2);            
        end
    end
end

%extract the max eigenvlaue of matrix A and C from the normalized eigenvextors
for i=1:group
    [eigvec_a(:,:,i),val_a(:,:,i)]=eig(A(:,:,i));
    [eigvec_c(:,:,i),val_c(:,:,i)]=eig(C(:,:,i));
    temp_a=A(:,:,i)*eigvec_a(:,:,i)/4;
    lamda_a(i)=max(temp_a(1,:)./eigvec_a(1,:,i));
    temp_c=C(:,:,i)*eigvec_c(:,:,i)/4;
    lamda_c(i)=max(temp_c(1,:)./eigvec_c(1,:,i));
end

%%generate related matrixs to simulate no huamn movements
clear i
for ii=1:N_tot
    if ii<=3
        H_ori_no(:,ii)=rand(N,1)+i*rand(N,1);  %generate 10 random channel measurement sampel
    else
        H_ori_no(:,ii)=H_ori_no(:,ii-3)+rand+i*rand;
    end
end

for i=1:group
    m=i;
    for j=1:N1
    H_no(:,j,i)=H_ori_no(:,m);  %put N1 samples in each time slot
    m=m+1;
    end
end

%calculate the correlation matrix A and C for amplitudes and phase of the CSI measurements
for ii=1:group
    for i=1:N
        for j=1:N
                tem=corrcoef(real(H_no(i,:,ii)),real(H_no(j,:,ii)));
                A_no(i,j,ii)=tem(1,2);
                tem=corrcoef(imag(H_no(i,:,ii)),imag(H_no(j,:,ii)));  
                C_no(i,j,ii)=tem(1,2);            
        end
    end
end

%extract the max eigenvlaue of matrix A and C from the normalized eigenvextors
for i=1:group
    [eigvec_a_no(:,:,i),val_a_no(:,:,i)]=eig(A_no(:,:,i));
    [eigvec_c_no(:,:,i),val_c_no(:,:,i)]=eig(C_no(:,:,i));
    tem_a=A_no(:,:,i)*eigvec_a_no(:,:,i)/4;
    lamda_a_no(i)=max(tem_a(1,:)./eigvec_a_no(1,:,i));
    tem_c=C_no(:,:,i)*eigvec_c_no(:,:,i)/4;
    lamda_c_no(i)=max(tem_c(1,:)./eigvec_c_no(1,:,i));
end


figure;
plot(lamda_a,lamda_c,'o');
xlabel('Normalized maximum eigenvalue \lambda_{A}');
ylabel('Normalized maximum eigenvalue \lambda_{C}');
grid on
hold on
plot(lamda_a_no,lamda_c_no,'*');
legend('movement','no movement');

figure;
plot(g,lamda_a);
hold on;
plot(g,lamda_a_no);
legend('movement \lambda_{A}','no movement \lambda_{A}');
xlabel('time');
ylabel('maximum eigenvalue');
grid on