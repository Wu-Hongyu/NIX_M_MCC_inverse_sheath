function [ v, position ] = pusher1D_leap_frog( v, position, q_m_ratio, E_particle ,simulation )
% 1D���������ǵ糡
% ��ֵ���ܵ͡����ڲ���������ʱʹ�ã�����Թ���ʱ
% ����ͨ�������޸ĵ�main��ʹ�ã���˱�۾�����main��һ�£�
% ���������޸�
        v(:,1)=v(:,1)+q_m_ratio*simulation.dt*E_particle;% ��ʱ��ڵ�
        position=position+v(:,1)*simulation.dt;% ��ʱ��ڵ�
end