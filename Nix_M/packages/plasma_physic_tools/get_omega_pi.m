function [ omega_pi ] = get_omega_pi( ni, Z, A)
% ����������徲����Ƶ��
% n: m^-3
% Z: �����
% A: ������
% omega_p: rad/s
constants=get_constants();
q=abs(Z)*constants.e;
omega_pi=sqrt(ni*q^2/(constants.eps0*A*constants.mH));
end

