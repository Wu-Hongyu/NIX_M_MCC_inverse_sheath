close all;clc;

load('��ԭ���ܶ�1e+16.mat');
figure(9)
plot(x_final,rhoH1n/-constants.e,'-b','LineWidth',2)
hold on
plot(x_final,rhoH/constants.e,'-r','LineWidth',2)
hold on
plot(x_final,rhoe/-constants.e,'-k','LineWidth',2)
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�����ܶȣ���ԭ���ܶ�1e+16��');
L1=legend('H-','H+','e');
set(L1,'location','north');
set(L1,'AutoUpdate','off');
figure(7)
plot(x_final,u_final,'-b','LineWidth',2);
hold on
figure(8)
plot(x_final,rho,'-b','LineWidth',2);
hold on
figure(14)
semilogy(x_final,rhoH1n/-constants.e,'-b','LineWidth',2);
hold on

load('��ԭ���ܶ�1e+17.mat');
figure(10)
plot(x_final,rhoH1n/-constants.e,'-b','LineWidth',2)
hold on
plot(x_final,rhoH/constants.e,'-r','LineWidth',2)
hold on
plot(x_final,rhoe/-constants.e,'-k','LineWidth',2)
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�����ܶȣ���ԭ���ܶ�1e+17��');
L1=legend('H-','H+','e');
set(L1,'location','north');
set(L1,'AutoUpdate','off');
figure(7)
plot(x_final,u_final,'-r','LineWidth',2);
hold on
figure(8)
plot(x_final,rho,'-r','LineWidth',2);
hold on
figure(14)
semilogy(x_final,rhoH1n/-constants.e,'-r','LineWidth',2);
hold on

load('��ԭ���ܶ�1e+18.mat');
figure(11)
plot(x_final,rhoH1n/-constants.e,'-b','LineWidth',2)
hold on
plot(x_final,rhoH/constants.e,'-r','LineWidth',2)
hold on
plot(x_final,rhoe/-constants.e,'-k','LineWidth',2)
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�����ܶȣ���ԭ���ܶ�1e+18��');
L1=legend('H-','H+','e');
set(L1,'location','north');
set(L1,'AutoUpdate','off');
figure(7)
plot(x_final,u_final,'-k','LineWidth',2);
hold on
figure(8)
plot(x_final,rho,'-k','LineWidth',2);
hold on
figure(14)
semilogy(x_final,rhoH1n/-constants.e,'-k','LineWidth',2);
hold on

load('��ԭ���ܶ�1e+19.mat');
figure(12)
plot(x_final,rhoH1n/-constants.e,'-b','LineWidth',2)
hold on
plot(x_final,rhoH/constants.e,'-r','LineWidth',2)
hold on
plot(x_final,rhoe/-constants.e,'-k','LineWidth',2)
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�����ܶȣ���ԭ���ܶ�1e+19��');
L1=legend('H-','H+','e');
set(L1,'location','north');
set(L1,'AutoUpdate','off');
figure(7)
plot(x_final,u_final,'-g','LineWidth',2);
hold on
figure(8)
plot(x_final,rho,'-g','LineWidth',2);
hold on
figure(14)
semilogy(x_final,rhoH1n/-constants.e,'-g','LineWidth',2);
hold on

load('��ԭ���ܶ�1e+20.mat');
figure(13)
plot(x_final,rhoH1n/-constants.e,'-b','LineWidth',2)
hold on
plot(x_final,rhoH/constants.e,'-r','LineWidth',2)
hold on
plot(x_final,rhoe/-constants.e,'-k','LineWidth',2)
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�����ܶȣ���ԭ���ܶ�1e+20��');
L1=legend('H-','H+','e');
set(L1,'location','north');
set(L1,'AutoUpdate','off');

figure(7)
plot(x_final,u_final,'-y','LineWidth',2);
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('\phi [V]');
L1=legend('��ԭ���ܶ�1e+16','��ԭ���ܶ�1e+17','��ԭ���ܶ�1e+18','��ԭ���ܶ�1e+19','��ԭ���ܶ�1e+20');
set(L1,'location','south');
set(L1,'AutoUpdate','off');

figure(8)
plot(x_final,rho,'-y','LineWidth',2);
hold on
axis([0,simulation.Lx*0.25,-inf,inf]);
xlabel('x [m]')
ylabel('������ܶ�\rho [C/m^3]');
L1=legend('��ԭ���ܶ�1e+16','��ԭ���ܶ�1e+17','��ԭ���ܶ�1e+18','��ԭ���ܶ�1e+19','��ԭ���ܶ�1e+20');
set(L1,'location','south');
set(L1,'AutoUpdate','off');

figure(14)
semilogy(x_final,rhoH1n/-constants.e,'-y','LineWidth',2);
hold on
axis([0,simulation.Lx*0.125,-inf,inf]);
xlabel('x [m]')
ylabel('�⸺�����ܶ�[/m^3]');
L1=legend('��ԭ���ܶ�1e+16','��ԭ���ܶ�1e+17','��ԭ���ܶ�1e+18','��ԭ���ܶ�1e+19','��ԭ���ܶ�1e+20');
set(L1,'location','north');
set(L1,'AutoUpdate','off');
