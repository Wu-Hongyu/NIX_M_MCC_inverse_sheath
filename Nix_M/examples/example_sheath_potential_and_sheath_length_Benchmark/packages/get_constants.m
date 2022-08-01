function [ constants ] = get_constants(  )
% ���� �洢���� �Ľṹ��
%--------������--------------------------------
constants.e=1.6022e-19;%������� [C]
constants.me=9.1094e-31;%�������� [kg]
constants.mH=1.6726e-27;%�������� [kg]
constants.c=2.9979e8;%��չ���[m/s]
constants.eps0=8.8542e-12;%��ս�糣�� [F/m]
constants.mu0=4*pi*1e-7;%��մŵ��� [H/m]
constants.kB=1.3807e-23;%������������ [ J/K]
constants.mH1n=1.6726e-27;%������������ [kg]
constants.mH2=3.3474e-27;%������������ [kg]
constants.mH3=5.0200e-27;

constants.q_m_ration_e=-constants.e/constants.me; %���Ӻ��ʱ� [C/kg]
constants.q_m_ration_H=constants.e/constants.mH; %���Ӻ��ʱ� [C/kg]
constants.q_m_ration_H1n=-constants.e/constants.mH1n; %�������Ӻ��ʱ� [C/kg]

constants.q_m_ratio_e=-constants.e/constants.me; %���Ӻ��ʱ� [C/kg]
constants.q_m_ratio_H=constants.e/constants.mH; %���Ӻ��ʱ� [C/kg]
constants.q_m_ratio_H1n=-constants.e/constants.mH1n; %�������Ӻ��ʱ� [C/kg]
constants.q_m_ratio_Hp=constants.e/constants.mH; 
constants.q_m_ratio_H2p=constants.e/constants.mH2; 
constants.q_m_ratio_H3p=constants.e/constants.mH3; 
constants.q_m_ratio_Hn=-constants.e/constants.mH; 
%--------������--------------------------------
end
