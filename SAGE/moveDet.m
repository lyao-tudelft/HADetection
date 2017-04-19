function [ cor_a, cor_c ] = moveDet( CSI, N1 )
% Compute parameters( maximum eigenvalues of correlation matrix) used to 
% detect movement.
%
% CSI: the channel state information
% N1: sliding window length for correlation
%
% cor_a: maximum eigenvalues of correlation matrix on amplitude
% cor_c: maximum eigenvalues of correlation matrix on phase

output_ori_no = CSI;
[~, nCSI] = size(CSI);
group = nCSI - N1;

for i=1:group
    m=i;
    for j=1:N1
    output_no(:,j,i)=output_ori_no(:,m);  %put N1 samples in each time slot
    m=m+1;
    end
end

%calculate the correlation matrix A and C for amplitudes and phase of the CSI measurements
for ii=1:group
    for i=1:N1
        for j=1:N1
                tem=corrcoef(real(output_no(:,i,ii)),real(output_no(:,j,ii)));
                A_no(i,j,ii)=tem(1,2);
                tem=corrcoef(imag(output_no(:,i,ii)),imag(output_no(:,j,ii)));  
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

cor_a = lamda_a_no;
cor_c = lamda_c_no;

end

