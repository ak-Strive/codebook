function [codebook_size] = calculate_codebook_size(P1,P2,Delta)


Xmax1=P1(1); Xmin1=P1(2); Ymax1=P1(3); Ymin1=P1(4); Zmax1=P1(5); Zmin1=P1(6);
Xmax2=P2(1); Xmin2=P2(2); Ymax2=P2(3); Ymin2=P2(4); Zmax2=P2(5); Zmin2=P2(6);
Xdelta1=Delta(1);Ydelta1=Delta(2);Zdelta1=Delta(3);Xdelta2=Delta(4);Ydelta2=Delta(5);Zdelta2=Delta(6);

% Xgrid1 = linspace(Xmin1,Xmax1,Xnum1); Ygrid1 = linspace(Ymin1,Ymax1,Ynum1); Zgrid1 = linspace(Zmin1,Zmax1,Znum1);
% Xgrid2 = linspace(Xmin2,Xmax2,Xnum2); Ygrid2 = linspace(Ymin2,Ymax2,Ynum2); Zgrid2 = linspace(Zmin2,Zmax2,Znum2);

Xgrid1=[Xmin1:Xdelta1:Xmax1]; Ygrid1=[Ymin1:Ydelta1:Ymax1]; Zgrid1=[Zmin1:Zdelta1:Zmax1]; 
Xgrid2=[Xmin2:Xdelta2:Xmax2]; Ygrid2=[Ymin2:Ydelta2:Ymax2]; Zgrid2=[Zmin2:Zdelta2:Zmax2]; 

Xnum1=length(Xgrid1); Ynum1=length(Ygrid1); Znum1=length(Zgrid1);
Xnum2=length(Xgrid2); Ynum2=length(Ygrid2); Znum2=length(Zgrid2);

codebook_size=Xnum1*Ynum1*Znum1*Xnum2*Ynum2*Znum2;
