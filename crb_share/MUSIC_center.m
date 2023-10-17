function [phi_est] = MUSIC_center(Y, Ns) %Ns:receiver antenna number
    Ry = Y*Y';
    [SR,~,~] = svd(Ry);
    UN = SR(:,2:Ns);

    J = 2*1e3;
    theta = -1:2/J:1;
    theta = theta(1:J);
    atm = ((0:(Ns-1))- (Ns-1)/2)';
    E = exp(1j*pi*atm*theta);
    P = sum(abs(E'*UN).^2,2);

    [~,ind] = sort(abs(P),'ascend');
    phi_est = theta(ind(1));
end

