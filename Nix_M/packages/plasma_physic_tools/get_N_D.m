function [ N_D ] = get_N_D( n, T )
% ����°�����������Ŀ�����ڴ���������ӣ�
% N_D: ��
% n: m^-3
% T: eV
constants=get_constants();
N_D=sqrt((constants.eps0*T/constants.e)^3/n)*4*pi/3;
% ע�⣬��SI��λ�ƣ�T in K���ı��ʽ��ͬ
end

