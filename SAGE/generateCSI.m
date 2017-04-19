function [ out ] = generateCSI(  N, N_tot, move )
% Generate CSIs of time-varying or stationary channel.
%
% N: number of path
% N_tot: number of CSIs
% move: time-varying if '1', stationary otherwise
%
% out: output CSI

total_N=1400;
%N_tot=size(theta,2);
frequency=1.4000e+10; %Ghz
delay=sort(randi(100,1,N));  %delay in times, unit is 'ns'
delay=floor((delay*10^-9)*frequency);   %delay at the correspongding index

n=1:N_tot;
H_ori=zeros(N,N_tot);
output_ori=rand(total_N,N_tot)/100+i*(rand(total_N,N_tot)/100); %add random noise to original output
output_ori_no=output_ori;
% output_ori = zeros(total_N,N_tot);
% output_ori_no = zeros(total_N,N_tot);
H_ori_no=zeros(N,length(n));
diff=zeros(10,10);
near=zeros(10,2);

if move == 1
    %%generate unrelated matrixs to simulate huamn movements
    for ii=1:N_tot
        if ii==1
            for j=1:10
                H_ori(j,1)=(2/j)*rand+i*rand;  %generate the first group of CSI of 10 paths
                output_ori(delay(j),1)=H_ori(j,1);
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
                if kk==0
                    R(4)=0;
                end
                R(5)=0;
                Rc(:,ii-1)=R';   %the number of path changed each time
                for j=1:10
                    if any(R==j)
                        H_ori(j,ii)=(2/j)*rand+i*rand;  %the paths change
                        difference=abs(H_ori(j,ii))-abs(H_ori(j,ii-1));
                        dis_delay=floor(delay(j)+difference*(total_N/N/2));
                        if dis_delay<=0
                            dis_delay=1;
                        end
                        if dis_delay>=1400
                            dis_delay=1400;
                        end
                        output_ori(dis_delay,ii)=H_ori(j,ii);
                    else
                        H_ori(j,ii)=H_ori(j,ii-1)+randn/10+i*(randn/10);  %the path
                        output_ori(delay(j),ii)=H_ori(j,ii);
                    end
                end
            end
        end
    end

    out = output_ori;
    
else
    %%generate related matrixs to simulate no huamn movements
    clear i
    for ii=1:N_tot
        if ii==1
            for j=1:10
            H_ori_no(j,1)=(2/j)*rand+i*rand; %generate the first group of CSI of 10 paths
            output_ori_no(delay(j),1)=H_ori_no(j,1);
            end
        else
            for j=1:10
            H_ori_no(j,ii)=H_ori_no(j,1)+randn/7+i*(randn/7);
            output_ori_no(delay(j),ii)=H_ori_no(j,ii);
            end
        end
    end
    out = output_ori_no;
    
end

end

