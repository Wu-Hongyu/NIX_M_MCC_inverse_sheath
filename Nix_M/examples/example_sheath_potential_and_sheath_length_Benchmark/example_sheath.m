%1D PIC 鞘层电压与长度benchmark

%% 初始化
clear
close all;

%文件路径
addpath('./packages')
%parameter
Te=2;
n0=4E17;
potential=zeros(size(Te));
length=zeros(size(n0));

for i=1:size(Te,2)
    [sheath_potential,sheath_length]=get_sheath(Te(i),n0(i),i);
    potential(i)=sheath_potential;
    length(i)=sheath_length;
end

figure;
hold on
h1=plot(Te,potential,'m','LineWidth',3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3D-PIC modelling of a low temperature plasma sheath with wall emission of negative particles and its application to NBI sources
%
the_1D_PIC_13_=[0.498737,	1.22222;...
    0.997475	2.51852;...
    1.99495	5.07407;...
    3.99621	10.037];

the_Onix=[            0.498737	1.35185;...
    0.991162	2.62963;...
    2.99242	7.33333;...
    4.99369	12.1111];

the_Analytic_37_=[0.0505051	0.037037;...
    4.72854	12.9815];

the_Analytic_39_=[0.0441919	0;...
    5.03157	12.9815];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2=plot(the_1D_PIC_13_(:,1),the_1D_PIC_13_(:,2),'--b','LineWidth',3);
h3=plot(the_Analytic_37_(:,1),the_Analytic_37_(:,2),'linestyle','-','Color',[255/256 128/256 64/256],'LineWidth',3);
h4=plot(the_Analytic_39_(:,1),the_Analytic_39_(:,2),'r','LineWidth',3);
h5=plot(the_Onix(:,1),the_Onix(:,2),'--g','LineWidth',3);
scatter(the_1D_PIC_13_(:,1),the_1D_PIC_13_(:,2),60,'b','LineWidth',3)
scatter(the_Onix(:,1),the_Onix(:,2),60,'g','LineWidth',3)
scatter(Te,potential,60,'m','LineWidth',3)
legend([h1,h2,h3,h4,h5],'Nix_M_','1D_PIC_13_','Analytic_37_','Analytic_39_','Onix','Location','NorthWest')
hold off

figure;

L1=loglog(n0,length,'m','LineWidth',3);
hold on
scatter(n0,length,'m','LineWidth',3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3D-PIC modelling of a low temperature plasma sheath with wall emission of negative particles and its application to NBI sources
%
the_Analytic_39_L=[8.197999e+013	0.00767283;...
    1.990824e+017	1.518504e-004];

the_Onix_L=[1.005638e+014	0.00651867;
    9.978768e+014	0.00206139;
    9.994654e+015	6.720996e-004;
    1.001125e+017	2.305854e-004;];

the_1D_PIC_L=[1.005735e+014	0.00700056;
    1.007267e+015	0.00216912;
    9.994516e+015	6.652867e-004;
    1.001002e+017	2.103821e-004;];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2=loglog(the_Analytic_39_L(:,1),the_Analytic_39_L(:,2),'r','LineWidth',3);

L4=plot(the_1D_PIC_L(:,1),the_1D_PIC_L(:,2),'--b','LineWidth',3);

L3=plot(the_Onix_L(:,1),the_Onix_L(:,2),'--g','LineWidth',3);

scatter(the_Onix_L(:,1),the_Onix_L(:,2),'g','LineWidth',3);
scatter(the_1D_PIC_L(:,1),the_1D_PIC_L(:,2),'b','LineWidth',3);

legend([L1,L2,L3,L4],'Nix_M_','Analytic_39_','Onix','1D_PIC_13_','Location','NorthEast')
hold off