clear

%load('theta.mat');  %results from sage saved in file theta.mat

N=10; %number of CSIs in each sample
N1=3; %number of samples in each time slot/window
N_tot=50;  %total number of channel measurement sample
%N_tot=size(theta,2);
group=N_tot-N1;

n=1:N_tot;
g=1:group;
H_ori=zeros(N,length(n));
H_ori_no=zeros(N,length(n));
H=zeros(N,N1,group);
A=zeros(N1,N1,group-1);
C=zeros(N1,N1,group-1);
A_no=zeros(N1,N1,group-1);
C_no=zeros(N1,N1,group-1);
lamda_a=zeros(1,group);
diff=zeros(10,10);
near=zeros(10,2);
Rc=zeros(5,group);

%%generate unrelated matrixs to simulate huamn movements
for ii=1:N_tot
    if ii==1
        for j=1:10
            H_ori(j,1)=(2/j)*rand+i*rand;  %generate the first group of CSI of 10 paths
        end
%if the difference of amplitudes of two received signals from two paths
%is less than 0.2, we reagrd these two path are so near that the change in
%one path wil have an impact on the other one
        k=1;
        for j=1:10
            for jj=1:10
                diff(j,jj)= abs(H_ori(j)-H_ori(jj));
                if j>=jj
                diff(j,jj)=0;
                end
                if diff(j,jj)>0 && diff(j,jj)<=0.2
                    near(k,:)=[j jj];
                    k=k+1;
                end
            end
        end
%with huamn activities, randomly chooseing three or four paths to change the rest keep unchanged
%maybe one more path will also be affected because of interference or
%Doppler effect
    else
        r=rand;
        if r>=0.5
            R=randperm(10,4); 
            kk=0;
                for i1=1:2
                    for j1=1:10
                        if kk==0
                            if any (R==near(j1,i1))
                                kk=1;
                                if i1==1
                                    R(5)=near(j1,2);
                                else
                                    R(5)=near(j1,1);  
                                end
                            end
                        end
                    end
                end
                if kk==0
                    R(5)=0;
                end
                Rc(:,ii-1)=R';
            for j=1:10
                if any(R==j)
                    H_ori(j,ii)=(2/j)*rand+i*rand;         
                else
                    H_ori(j,ii)=H_ori(j,ii-1)+abs(randn/10);
                end
            end
        else
            R=randperm(10,3);
            kk=0;
            for i1=1:2
                for j1=1:10
                    if kk==0
                        if any (R==near(j1,i1))
                            kk=1;
                            if i1==1
                                R(4)=near(j1,2);
                            else
                                R(4)=near(j1,1);   %to find if there is one more path being affected or not, due to Doppler
                            end
                        end
                    end
                end
            end
            R(5)=0;
            Rc(:,ii-1)=R';   %the number of path changed each time
            for j=1:10
                if any(R==j)
                    H_ori(j,ii)=(2/j)*rand+i*rand;  %the paths change
                else
                    H_ori(j,ii)=H_ori(j,ii-1)+randn/10;  %the path 
                end
            end
        end
    end
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
    for i=1:N1
        for j=1:N1
                temp=corrcoef(real(H(:,i,ii)),real(H(:,j,ii)));
                A(i,j,ii)=temp(1,2);
                temp=corrcoef(imag(H(:,i,ii)),imag(H(:,j,ii)));  
                C(i,j,ii)=temp(1,2);            
        end
    end
end

%extract the max eigenvlaue of matrix A and C from the normalized eigenvextors
for i=1:group
    [eigvec_a(:,:,i),val_a(:,:,i)]=eig(A(:,:,i));
    [eigvec_c(:,:,i),val_c(:,:,i)]=eig(C(:,:,i));
    temp_a=A(:,:,i)*eigvec_a(:,:,i)/N1;
    lamda_a(i)=max(temp_a(1,:)./eigvec_a(1,:,i));
    temp_c=C(:,:,i)*eigvec_c(:,:,i)/N1;
    lamda_c(i)=max(temp_c(1,:)./eigvec_c(1,:,i));
end

%%generate related matrixs to simulate no huamn movements
clear i
for ii=1:N_tot
    if ii==1
        for j=1:10
        H_ori_no(j,1)=(2/j)*rand+i*rand; %generate the first group of CSI of 10 paths
        end
    else
        for j=1:10
        H_ori_no(j,ii)=H_ori_no(j,1)+randn/10+i*(randn/10);
        end
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
    for i=1:N1
        for j=1:N1
                tem=corrcoef(real(H_no(:,i,ii)),real(H_no(:,j,ii)));
                A_no(i,j,ii)=tem(1,2);
                tem=corrcoef(imag(H_no(:,i,ii)),imag(H_no(:,j,ii)));  
                C_no(i,j,ii)=tem(1,2);            
        end
    end
end

%extract the max eigenvlaue of matrix A and C from the normalized eigenvextors
for i=1:group
    [eigvec_a_no(:,:,i),val_a_no(:,:,i)]=eig(A_no(:,:,i));
    [eigvec_c_no(:,:,i),val_c_no(:,:,i)]=eig(C_no(:,:,i));
    tem_a=A_no(:,:,i)*eigvec_a_no(:,:,i)/N1;
    lamda_a_no(i)=max(tem_a(1,:)./eigvec_a_no(1,:,i));
    tem_c=C_no(:,:,i)*eigvec_c_no(:,:,i)/N1;
    lamda_c_no(i)=max(tem_c(1,:)./eigvec_c_no(1,:,i));
end


figure;
plot(lamda_a,lamda_c,'o');
xlabel('Normalized maximum eigenvalue \lambda_{A}');
ylabel('Normalized maximum eigenvalue \lambda_{C}');
title('Normalized maximum eigenvalues of amplitude and phase correlation matrices');
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
title('Normalized maximum eigenvalues of amplitude');
grid on

figure;
plot(g,lamda_c);
hold on;
plot(g,lamda_c_no);
legend('movement \lambda_{A}','no movement \lambda_{C}');
xlabel('time');
ylabel('maximum eigenvalue');
title('Normalized maximum eigenvalues of phase');
grid on
