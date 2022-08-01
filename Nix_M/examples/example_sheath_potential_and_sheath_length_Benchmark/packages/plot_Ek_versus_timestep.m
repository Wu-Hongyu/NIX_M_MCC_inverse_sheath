function [ ] = plot_Ek_versus_timestep( vmean1,vmean2,q_m_ration1, q_m_ration2, name1,name2 ,dptime )
% ��ͼ��ʾ ���Ӷ���-ʱ�䲽
Ek=@(v,q_m_ratio) v.*v*3/(2*abs(q_m_ratio)); % �����Ժ��ᵽ����ȥ
plot(Ek(vmean1,q_m_ration1),'-b')
hold on
plot(Ek(vmean2,q_m_ration2),'--r')
xlabel(['ʱ�� [' num2str(dptime) 'ʱ��]']);
ylabel('3*x��Ek [eV]');
% title('���Ӷ���ʱ���ݻ�');
L1=legend(name1,name2);
set(L1,'location','best');
hold off
end