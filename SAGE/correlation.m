% Generate CSI in stationary and time-varying channels, compute parameters
% (maximum eigenvalues) to detect movement
clear

%load('theta.mat');  %results from sage saved in file theta.mat

%% Setting parameters
fs = 3.2e6;
N=10; %number of path
total_N=2000;
N1=3; %number of samples in each time slot/window
N_tot=50;  %total number of channel measurement sample
%N_tot=size(theta,2);
group=N_tot-N1;
frequency=fs; %Ghz
delay=sort(rand(1,N))*total_N/fs;  %delay in times
delay=floor(delay*frequency);   %delay at the correspongding index

n=1:N_tot;
g=1:group;
H_ori=zeros(N,N_tot);
H_ori_no=zeros(N,length(n));
A=zeros(N1,N1,group-1);
C=zeros(N1,N1,group-1);
A_no=zeros(N1,N1,group-1);
C_no=zeros(N1,N1,group-1);
lamda_a=zeros(1,group);
diff=zeros(10,10);
near=zeros(10,2);
Rc=zeros(5,group);

%% Movement Detection
%%generate unrelated matrixs to simulate huamn movements
output_ori = generateCSI(N, N_tot, 1);
amp = cell(N_tot,1);
[ lamda_a, lamda_c ] = moveDet( output_ori, N1 );

%%generate related matrixs to simulate no huamn movements
output_ori_no = generateCSI(N, N_tot, 0);
[ lamda_a_no, lamda_c_no ] = moveDet( output_ori_no, N1 );

%% Plot results
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