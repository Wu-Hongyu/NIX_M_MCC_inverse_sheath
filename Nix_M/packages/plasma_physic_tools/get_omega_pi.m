function [ omega_pi ] = get_omega_pi( ni, Z, A)
% 计算等离子体静电振荡频率
% n: m^-3
% Z: 电荷数
% A: 质量数
% omega_p: rad/s
constants=get_constants();
q=abs(Z)*constants.e;
omega_pi=sqrt(ni*q^2/(constants.eps0*A*constants.mH));
end

