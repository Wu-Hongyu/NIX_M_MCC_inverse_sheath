function [ debye_length ] = get_debye_length( n,T )
% ����°ݳ��ȣ����ڴ���������ӣ�
% n: m^-3
% T: eV
% debye_length: m
constants=get_constants();
debye_length=sqrt(constants.eps0*T/(constants.e*n));
% �����ڵķ�ʽ�����ӷ�ĸ��ȥ��һ��constants.e
end