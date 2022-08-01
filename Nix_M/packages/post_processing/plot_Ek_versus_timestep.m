function [ ] = plot_Ek_versus_timestep( vmean1,vmean2,vmean3,vmean4,q_m_ration1, q_m_ration2, name1,name2 ,name3,name4,dptime )
% 画图显示 粒子动能-时间步
Ek=@(v,q_m_ratio) v.*v/(2*abs(q_m_ratio)); % 可能以后提到外面去
plot(Ek(vmean1,q_m_ration1),'-b')
hold on
plot(Ek(vmean2,q_m_ration2),'--r')
plot(Ek(vmean3,q_m_ration2/2),'--g')
plot(Ek(vmean4,q_m_ration2/3),'--m')


xlabel(['时间 [' num2str(dptime) '时步]']);
ylabel('3*x向Ek [eV]');
% title('粒子动能时间演化');
L1=legend(name1,name2,name3,name4);
set(L1,'location','best');
hold off
end