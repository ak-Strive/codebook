clc; 
clear all
close all

M = 1; % the number of the antennas at the BS
N1 = 128; 
N2 = 4;
N = N1*N2; % the number of RIS elements
d = 0.5; % the antenna spacing lambda/2

A = 4;
delta = 0.25;
sample = 10;

SNR_dB = [-5:2:5];  
SNR_linear = 10.^(SNR_dB/10.);

% generate the far-field codebook 
UN1 = exp(1i*2*pi*[0:N1-1]'*d*[0:N1-1]*(2/N1));
UN2 = exp(1i*2*pi*[0:N2-1]'*d*[0:N2-1]*(2/N2));
far_codebook = kron(UN1,UN2);

% generate the near-field codebook
P1 = [1200*d,-1200*d,200*d,10*d,400*d,-400*d];
P2 = [1200*d,-1200*d,200*d,10*d,400*d,-400*d];
Delta = [100*d,100*d,100*d,100*d,100*d,100*d];
[near_codebook,~] = generate_near_field_codebook(N1,N2,d,P1,P2,Delta);

% generate the near-field hierachical codebook
Delta1 = Delta*A;
[near_codebook1,record] = generate_near_field_codebook(N1,N2,d,P1,P2,Delta1);

rate_far = zeros(sample,length(SNR_dB));
rate_near = zeros(sample,length(SNR_dB));
rate_near2 = zeros(sample,length(SNR_dB));
rate_opt = zeros(sample,length(SNR_dB));


for t = 1:sample
    t
% generate the channel from the BS to the RIS
[G,px1,py1,pz1] = generate_G_near_field_channel(N1,N2,P1);

% generate the channel from the RIS to the user
[hK,px2,py2,pz2] = generate_hr_near_field_channel(N1,N2,1,P2);

for s = 1:length(SNR_dB)
    s
    SNR = SNR_linear(s);
    Hc = diag(hK)*G; 
    
    %% far-field beam training 
    array_gain = 0;
    for i =1:length(far_codebook)
        array_gain=max(array_gain,abs(far_codebook(i,:)*Hc)^2);
    end
    rate_far(t,s) = log2(1 + SNR * array_gain);
    
    %% near-field beam training 
    array_gain = 0;
    for i =1:size(near_codebook,1)
        array_gain=max(array_gain,abs(near_codebook(i,:)*Hc)^2);
    end
    rate_near(t,s) = log2(1 + SNR * array_gain);
    
    %% hierarchical near-field beam training 
    array_gain = 0;
    max_index=-1;
    for i =1:size(near_codebook1,1)
        if abs(near_codebook1(i,:)*Hc)^2>array_gain
            max_index=i;
            array_gain=abs(near_codebook1(i,:)*Hc)^2;
        end
    end
    % generate the second-level codes
    P21=[record(max_index,1)+Delta1(1)/2,record(max_index,1)-Delta1(1)/2,record(max_index,2)+Delta1(2)/2,record(max_index,2)-Delta1(2)/2,record(max_index,3)+Delta1(3)/2,record(max_index,3)-Delta1(3)/2];
    P22=[record(max_index,4)+Delta1(4)/2,record(max_index,4)-Delta1(4)/2,record(max_index,5)+Delta1(5)/2,record(max_index,5)-Delta1(5)/2,record(max_index,6)+Delta1(6)/2,record(max_index,6)-Delta1(6)/2];
    
    near_codebook2 = generate_near_field_codebook(N1,N2,d,P21,P22,Delta1*delta);
    
    for i =1:size(near_codebook2,1)
        if abs(near_codebook2(i,:)*Hc)^2>array_gain
            array_gain=abs(near_codebook2(i,:)*Hc)^2;
        end
    end
    rate_near2(t,s) = log2(1 + SNR * array_gain);
    
    %% perfect CSI based beamforming
    wc_opt = exp(1j*phase(Hc'));
    array_gain = abs(wc_opt*Hc)^2;
    rate_opt(t,s) = log2(1 + SNR * array_gain);
end
end

figure('color',[1,1,1]); hold on; box on; grid on;
ha=gca;
plot(SNR_dB,mean(rate_far),'b^-', 'Linewidth', 1.6)
plot(SNR_dB,mean(rate_near),'ms-','Linewidth',1.6)
plot(SNR_dB,mean(rate_near2),'rd-','Linewidth',1.6)
plot(SNR_dB,mean(rate_opt),'k--','Linewidth', 1.6)
legend('Far-field beam training [17]','Proposed near-field beam training','Proposed hierarchical near-field beam training','Perfect CSI based beamforming')

xlabel('SNR (dB)')
ylabel('Achievable Rate (bit/s/Hz)')

% save results.mat rate_*

set(ha,'xtick',[-5:2:5])
set(ha,'ytick',[11:2:19])
axis([-5,5,11,19])

