function [ omega_pe ] = get_omega_pe(ne)
% ������ӵ�������Ƶ��
% n: m^-3
% omega_pe: rad/s
constants=get_constants();
omega_pe=sqrt(ne*constants.e*constants.e/(constants.eps0*constants.me));
end

