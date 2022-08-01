function [ ] = plot_vineV_versus_position( v1, v2, position1,position2, q_m_ration1, q_m_ration2, name1,name2 ,simulation )
% 画图显示 相空间
vineV=@(v,q_m_ratio) sign(v).*v.*v*3/(2*abs(q_m_ratio));
scatter(position1,vineV(v1(:,1),q_m_ration1),1,'b')
axis([0 simulation.Lx -inf inf]);
hold on
scatter(position2,vineV(v2(:,1),q_m_ration2),1,'r')
%         title('粒子相空间分布');
xlabel('x [m]');
ylabel('directed Ek [eV]');
L1=legend(name1,name2);
set(L1,'location','northwest');
hold off
end