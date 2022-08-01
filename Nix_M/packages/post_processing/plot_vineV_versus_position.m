function [ ] = plot_vineV_versus_position( v1, v2,v3,v4, position1,position2,position3,position4, q_m_ration1, q_m_ration2, name1,name2,name3,name4 ,simulation )
% 画图显示 相空间
yyaxis left
vineV=@(v,q_m_ratio) v.*v/(2*abs(q_m_ratio));
scatter(position1,vineV(v1,q_m_ration1),1,'b')
axis([0 simulation.Lx -inf inf]);
%hold on
yyaxis right
scatter(position2,vineV(v2,q_m_ration2),1,'r')
hold on
scatter(position3,vineV(v3,q_m_ration2/2),1,'g')
hold on
scatter(position4,vineV(v4,q_m_ration2/3),1,'m')
hold off
hold off
%         title('粒子相空间分布');
xlabel('x [m]');
ylabel('directed Ek [eV]');
L1=legend(name1,name2,name3,name4);
set(L1,'location','northwest');
hold off
end