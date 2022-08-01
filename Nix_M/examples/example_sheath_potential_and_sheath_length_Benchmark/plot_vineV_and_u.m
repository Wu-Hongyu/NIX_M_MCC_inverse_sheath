function [ ] = plot_vineV_and_u( vH_grid, u,  q_m_ration, constants ,Te)
% 画图显示 相空间
vineV=@(v,q_m_ratio)v.*v*3/(2*abs(q_m_ratio));
vH_grid_eV=vineV(vH_grid,q_m_ration);
yyaxis left
plot(vH_grid_eV,'-b','LineWidth',3)
axis([0 size(vH_grid,1) -inf inf]);
hold on
uB0=sqrt(constants.e*Te/constants.mH);
uB=vineV(uB0,q_m_ration);
line([0,size(vH_grid,1)],[uB,uB],'linestyle','-.','color','b','LineWidth',2)
xlabel('网格编号');
ylabel('H Ek [eV]');
hold off
yyaxis right
hold off
plot(u,'-r','LineWidth',2)
ylabel('u [V]');
% uB0=sqrt(constants.e*Te/constants.mH);
% uB=vineV(uB0,q_m_ration);
% line([0,size(vH_grid,1)],[uB,uB],'linestyle','-.','color','b','LineWidth',2)
L1=legend('离子动能','玻姆速度','电位');
set(L1,'location','northwest');
hold off
end