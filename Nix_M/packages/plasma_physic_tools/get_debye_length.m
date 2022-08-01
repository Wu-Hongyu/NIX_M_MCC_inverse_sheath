function [ debye_length ] = get_debye_length( n,T )
% 计算德拜长度（限于带单电荷粒子）
% n: m^-3
% T: eV
% debye_length: m
constants=get_constants();
debye_length=sqrt(constants.eps0*T/(constants.e*n));
% 根号内的分式，分子分母消去了一个constants.e
end