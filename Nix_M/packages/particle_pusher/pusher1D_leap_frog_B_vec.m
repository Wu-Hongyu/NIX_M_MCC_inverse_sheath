function [ v, position ] = pusher1D_leap_frog_B_vec( v, position, q_m_ratio, E_particle,B_particle,simulation )
% 1D���������ǵ糡��ʸ���������Ǵų�
% ��ֵ���ܵ͡����ڲ���������ʱʹ�ã�����Թ���ʱ
% ����ͨ�������޸ĵ�main��ʹ�ã���˱�۾�����main��һ�£�
% ���������޸�
        v(:,1)=v(:,1)+q_m_ratio*simulation.dt*(E_particle-v(:,2).*B_particle);% ��ʱ��ڵ�
        v(:,2)=v(:,2)+q_m_ratio*simulation.dt*v(:,1).*B_particle; %Bz
%         vpx=vpx+qm*(mat*Eg+vpy*Bz0-vpz*By0)*dt;
%         vpy=vpy+qm*(vpz*Bx0-vpx*Bz0)*dt;
%         vpz=vpz+qm*(vpx*By0-vpy*Bx0)*dt;
        position=position+v(:,1)*simulation.dt;% ��ʱ��ڵ�
end