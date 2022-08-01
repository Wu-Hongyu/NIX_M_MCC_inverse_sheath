function [ N_D ] = get_N_D( n, T )
% 计算德拜球内粒子数目（限于带单电荷粒子）
% N_D: 个
% n: m^-3
% T: eV
constants=get_constants();
N_D=sqrt((constants.eps0*T/constants.e)^3/n)*4*pi/3;
% 注意，与SI单位制（T in K）的表达式不同
end

