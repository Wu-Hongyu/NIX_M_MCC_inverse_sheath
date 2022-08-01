function [k]=get_k_e_heavy(xsec, Te)
% 无电磁场时，Maxwellian分布电子撞击重粒子反应的速率系数
constants=get_constants();
%截面对Maxwellian分布积分
get_k=@(E,a,Te) trapz(E,a.*E.*exp(-E/Te))/Te^1.5*...
    sqrt(8*abs(constants.q_m_ratio_e)/pi); 
k=get_k(xsec(:,1),xsec(:,2),Te);
end

%PengChen: test ok 2020/12/24 17:13:44 
% 通过get_k_goal_Te_range.m进行测试