clear;
% 运行三次，theat_IT = 80, 40, 0 degree
tic
Nb = 64; % BS antenna number
M = 64; % RIS reflecting number
Ms = 8; % RIS sensing number
fc = 28*1e9; % carrier frequencyd
lambda = 3*1e8/fc;
d_br = 30;
d_ru = 10;
d_rt = 5;
alpha_br2 = (lambda/4/pi/d_br)^2; %alpha_g
alpha_ru2 = (lambda/4/pi/d_ru)^2; %alpha_h
kappa = 7; %dBsm
kappa = 10^(kappa/10);
alpha_rt2 = lambda^2*kappa/64/pi^3/d_rt^4; %alpha_s
%%
theta_BI= -30; 
theta_IT = 40; % target angle
theta_IT_rad = theta_IT*pi/180;
a_s = exp(1j*pi*((0:Ms-1)-(Ms-1)/2)'*sind(theta_IT));
q_r = exp(1j*pi*((0:M-1)- (M-1)/2)'* (sind(theta_BI)- sind(theta_IT)));
%%
T = 1000; % 相干时间 
sigma2 = 1e-15;     % 噪声功率，可转化为-120dBm
ite_num = 100 ;     % independent realizations迭代数
method = 2;  % method 1: music estimation; method 2: MLE estimation
Pt_list_dB = 0:2:30; %dBm
Pt_list = 10.^(Pt_list_dB/10-3);
pt_len = length(Pt_list);
%% crb_theory
steering_ind = -(M +1-2*(1:M )')/2;
% q_steering = exp(1j*pi*(sind(theta_BI)- sind(theta_IT))*steering_ind);  
dao_q_steering = -1j*pi*cosd(theta_IT).*steering_ind.*q_r;
norm_q = (norm(dao_q_steering))^2;

steering_ind = -(Ms+1-2*(1:Ms)')/2;
b_steering = exp(1j*pi*sind(theta_IT)*steering_ind);
dao_b_steering = 1j*pi*cosd(theta_IT).*steering_ind.*b_steering;
norm_b = (norm(dao_b_steering))^2;
%%
NMSE = zeros(ite_num, pt_len);
for i=1:pt_len
    fprintf('Progress:%.2f%%\n',i/pt_len*100);
    Pt = Pt_list(i);
    tau = M; % time of beam sweeping, codebook size
    %% codebook
    eta = -1+(2*(1:tau)-1)/tau; 
    D = exp(1j*pi*((0:M-1)- (M-1)/2)'*eta);
    %%
    for ite = 1:ite_num
        noise = sqrt(sigma2/2)*(randn(Ms,tau)+1j*randn(Ms,tau)); 
        ytm = sqrt(Nb)*sqrt(Pt)*sqrt(alpha_br2)*sqrt(alpha_rt2)*a_s*transpose(q_r)*D;
        Y = ytm+noise;
        %% MUSIC
        if method == 1
            theta_est = MUSIC_center(Y, Ms);
    %         theta_est = MUSIC(Y, Ms);
            NMSE(ite,i) = (theta_IT*pi/180-asind(theta_est)*pi/180)^2;
        elseif method == 2
        %% MLE 使用极大似然估计
            J = 2*1e3;
            theta_list = -1:2/J:1;
            q_ind = -(M +1-2*(1:M )')/2;
            as_ind = -(Ms+1-2*(1:Ms)')/2;
            re_mle = zeros(J,1);
            for ind_j=1:J
                theta = theta_list(ind_j);
                as_tm = exp(1j*pi*as_ind*theta);
                q_tm = exp(1j*pi*q_ind*(sind(theta_BI)-theta));
                re_tm = as_tm'*Y*D'*conj(q_tm);
                re_mle(ind_j,1) = (abs(re_tm))^2;
            end
            [~,re_ind]=sort(re_mle,'descend');
            theta_est = asind(theta_list(re_ind(1)));
            NMSE(ite,i) = (theta_IT*pi/180-theta_est*pi/180)^2;
        end
    end
    CRB(i,1) = sigma2/2/tau/alpha_rt2/alpha_br2/ (Ms*norm_q+ M*norm_b)/Pt/Nb;
end
toc
a = sigma2/alpha_br2/alpha_rt2;
RMSE = sqrt(mean(NMSE));
CRB = sqrt(CRB);

figure; semilogy(Pt_list_dB,RMSE,'-b*','LineWidth',1.5);grid on;
hold on;semilogy(Pt_list_dB,CRB,'--b','LineWidth',1.5);
xlabel('Transmit Power (dBm)','FontSize',16);ylabel('RMSE','FontSize',16);

hold on;semilogy(Pt_list_dB,RMSE,'-mo','LineWidth',1.5); 
hold on;semilogy(Pt_list_dB,CRB,'--m','LineWidth',1.5);

hold on;semilogy(Pt_list_dB,RMSE,'-k+','LineWidth',1.5); 
hold on;semilogy(Pt_list_dB,CRB,'--k','LineWidth',1.5);

legend({'RMSE, \theta = 80^\circ', 'Root-CRB, \theta = 80^\circ',...
    'RMSE, \theta = 40^\circ', 'Root-CRB, \theta = 40^\circ',...
    'RMSE, \theta = 0^\circ', 'Root-CRB, \theta = 0^\circ'},'FontSize',13)
