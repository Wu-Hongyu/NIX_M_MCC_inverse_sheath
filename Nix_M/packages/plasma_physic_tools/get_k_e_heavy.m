function [k]=get_k_e_heavy(xsec, Te)
% �޵�ų�ʱ��Maxwellian�ֲ�����ײ�������ӷ�Ӧ������ϵ��
constants=get_constants();
%�����Maxwellian�ֲ�����
get_k=@(E,a,Te) trapz(E,a.*E.*exp(-E/Te))/Te^1.5*...
    sqrt(8*abs(constants.q_m_ratio_e)/pi); 
k=get_k(xsec(:,1),xsec(:,2),Te);
end

%PengChen: test ok 2020/12/24 17:13:44 
% ͨ��get_k_goal_Te_range.m���в���