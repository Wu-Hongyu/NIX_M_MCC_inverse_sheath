% -*- coding: utf-8 -*-
% ----------------------------------------------
%{
 * brief #Abstract  2018zhangģ�������������
 * Created 21:14:46 2022/07/14
 * author PengChen, HUST(peng_chen2016@hust.edu.cn)
 *
 * note #Detail

base·����f:\Project\2018zhang_H_inserve
�ĵ���20220714����-2018zhang���ʲ�ģ��.docx
���룺inverse - ����\code\

 * #TODO
%}
% ----------------------------------------------
addpath('../packages')
constants=get_constants();

%% ��ɢ����
dt=5.d-12; % ʱ�䲽��

ZLength=0.002; % �ܳ���
NxMax=513; % ���������
dx=ZLength/(NxMax-1);

InitDensity=4.0e17;
ParticlePerGrid=100;
Weight=InitDensity/ParticlePerGrid;

PB_VFactor=dx/dt;

%% ���Ӳ���
ne=4e17;
nHp=4e17;
nHn=3.52e17;
nH=1e20;

Te=2; % eV
THp=0.8; % eV
THn=0.8; % eV
TH=0.8; % eV
% Te*constants.e/constants.kB % 23210 K

% OES��ϵõ���ԭ���¶�
% 0.8eV	  2006Fantz
% 0.19eV:2.5eV=1:1	 2018Fantz

% ������ܶ�		D-E	��1e19m^-3	��ѹ����	2017Heinemann
% ��ԭ�ӱ���		
% E	0.2	OES	2006Fantz
% BUG D-E	0.3��0.1	OES	2018Fantz

% 2014Wimmer���������
%% ��ɽ�����ײ
THp_cold=0.01; % eV
% THp_cold*constants.e/constants.kB % 116.0426 K
% (0.5*constants.mH*PB_VFactor^2)/constants.e % 3.1859e+03 K���������ɶ���3�� �� 1000
% CXTemperature=0.11605; % ����Ϊ THp_cold/1000

% ʹ�ö����ٶ� ���������ɶ�
vHp=sqrt((THp*(3./2)*constants.e)/(0.5*constants.mH));
% VFactor=1./PB_VFactor;
% V =V/VFactor;
V=vHp;
Q=(-0.8221*log(V)+15.1262)*1.0e-20;
% f=V*Q*1.0e21*5; % ��֪2018Zhangʹ��1000f����1000*nH*V*Q
f=1*V*Q*nH;
P0=1-exp(-f*dt);
P1000=1-exp(-1000*f*dt);
P5=1-exp(-5*f*dt);

% 2018Zhang�����ͼ
Q_fun=@(V) (-0.8221*log(V)+15.1262)*1.0e-20;
x_T=THp/10:THp/10:10*THp;
% ʹ�ö����ٶ� ���������ɶ�
x_v=sqrt((x_T*(3./2)*constants.e)/(0.5*constants.mH));
y=Q_fun(x_v)*5;
figure
plot(x_T,y)

grid on
hold on
hev=[1.21E-01;1.66E-01;2.27E-01;3.12E-01;4.28E-01;5.87E-01;8.04E-01;1.10E+00;1.51E+00;2.07E+00;2.85E+00;
3.90E+00;5.35E+00;7.34E+00;1.01E+01;1.38E+01;1.89E+01;2.59E+01;3.56E+01;4.88E+01;6.69E+01;9.18E+01];
xec=[5.04E-19;4.65E-19;4.50E-19;4.46E-19;4.43E-19;4.38E-19;4.29E-19;4.17E-19;4.01E-19;3.84E-19;3.67E-19;
3.51E-19;3.36E-19;3.24E-19;3.13E-19;3.03E-19;2.95E-19;2.87E-19;2.80E-19;2.71E-19;2.63E-19;2.53E-19];
hev2=[1.00E-01;1.26E-01;1.58E-01;2.00E-01;2.51E-01;3.16E-01;3.98E-01;5.01E-01;6.31E-01;7.94E-01;1.00E+00;
1.26E+00;1.58E+00;2.00E+00;2.51E+00;3.16E+00;3.98E+00;5.01E+00;6.31E+00;7.94E+00;1.00E+01];
xec2=[2.19E-18;2.36E-18;2.05E-18;2.21E-18;2.06E-18;1.89E-18;1.82E-18;1.90E-18;1.65E-18;1.59E-18;1.63E-18;
1.51E-18;1.49E-18;1.45E-18;1.36E-18;1.36E-18;1.27E-18;1.26E-18;1.19E-18;1.17E-18;1.10E-18];
eV=0.5*constants.mH*V*V/constants.e;
Q2=interp1(hev,xec,eV);
f=1*V*Q2*nH;
% dt=1.7e-11;
Pp=1-exp(-f*dt);


plot(hev,xec)
ylabel('{\it\sigma}_{CEX}')
xlabel('{\itv}_{H^+} [eV]')
axis([0,10,-inf,inf])


%% �������
% % ʹ�����ٶ�Ϊ�����ٶ� ���������ɶ�
% vH=sqrt((TH*constants.e)/(0.5*constants.mH));
% vHp=sqrt((THp*constants.e)/(0.5*constants.mH));
% ʹ�ö����ٶ� ���������ɶ�
vH=sqrt((TH*(3./2)*constants.e)/(0.5*constants.mH));
vHp=sqrt((THp*(3./2)*constants.e)/(0.5*constants.mH));
Gamma_H=nH*vH/4;
Gamma_Hp=nHp*vHp/4;

% ��֪��ƽ�����ӵ��¶�
fY_T=@(T,RN_eta0,Eth_div_RE) RN_eta0*exp(-Eth_div_RE./T); % T in eV
% ����Cs�㣨����Լ1.5eV������
RN_eta0_Hion=0.3;
Eth_div_RE_Hion=2;
RN_eta0_Hatom=0.42;
Eth_div_RE_Hatom=1.05;
% ������
fY_T_Hatom=@(T) fY_T(T,RN_eta0_Hatom,Eth_div_RE_Hatom);
fY_T_Hion=@(T) fY_T(T,RN_eta0_Hion,Eth_div_RE_Hion);

% ��֪��ƽ���¶�
Gamma_Hn=Gamma_H*fY_T_Hatom(TH)+Gamma_Hp*fY_T_Hion(THp);
% 2.18e21


% %% 2014Wimmer����
% nH2=2e19;
% H_ratio=[0.2,0.4];
% TH=0.8;
% nH=nH2*H_ratio;
% 
% nion=[0.5,5]*1e17;
% Hp_ratio=1/2.5;
% nHp=nion*Hp_ratio;
% THp=0.8;
% 
% vH=sqrt((TH*constants.e)/(0.5*constants.mH));
% vHp=sqrt((THp*constants.e)/(0.5*constants.mH));
% Gamma_H=nH*vH/4;  % 1.0e+22 * [1.2380    2.4760]
% Gamma_Hp=nHp*vHp/4; % 1.0e+20 * [0.6190    6.1900]
% 
% % H-����from���ܶ�
% 
% % H-����from�������
% Gamma_Hn=Gamma_H*fY_T_Hatom(TH)+Gamma_Hp*fY_T_Hion(THp);

p=0.3;
T=1200;
n=p/(T*constants.kB)
