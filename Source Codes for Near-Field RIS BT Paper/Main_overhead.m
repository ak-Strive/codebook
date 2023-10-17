clc; 
clear all;
close all;

d=0.5;
A=2;
delta=0.25;

% generate the near-field codebook
P1=[1200*d,-1200*d,200*d,10*d,400*d,-400*d];
P2=[1200*d,-1200*d,200*d,10*d,400*d,-400*d];


xd_all=[80:10:150];

codebook_size = zeros(length(xd_all),1);
hierarchical_codebook_size = zeros(length(xd_all),1);

for i=1:length(xd_all)
    xd=xd_all(i);
    Delta=[xd*d,xd*d,xd*d,xd*d,xd*d,xd*d];
    P12=[-xd*d/2,xd*d/2,-xd*d/2,xd*d/2,-xd*d/2,xd*d/2];
    P22=[-xd*d/2,xd*d/2,-xd*d/2,xd*d/2,-xd*d/2,xd*d/2];
    codebook_size(i)=calculate_codebook_size(P1,P2,Delta);
    hierarchical_codebook_size(i)=calculate_codebook_size(P1,P2,Delta*A)+calculate_codebook_size(P12,P22,Delta*A*delta);
end

figure('color',[1,1,1]); hold on; box on; grid on;
ha=gca;
plot(xd_all,codebook_size,'ms-','Linewidth',1.6)
plot(xd_all,hierarchical_codebook_size,'rd-','Linewidth',1.6)
legend('Proposed near-field beam training','Proposed hierarchical near-field beam training')

xlabel('Sampling Step')
ylabel('Beamtraining Training Overhead')

