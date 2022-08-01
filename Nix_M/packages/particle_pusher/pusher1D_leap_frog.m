function [ v, position ] = pusher1D_leap_frog( v, position, q_m_ratio, E_particle ,simulation )
% 1D蛙跳法考虑电场
% 传值性能低。仅在不在意性能时使用，如测试功能时
% 可能通过复制修改到main中使用，因此表观尽量与main中一致，
% 减少所需修改
        v(:,1)=v(:,1)+q_m_ratio*simulation.dt*E_particle;% 半时间节点
        position=position+v(:,1)*simulation.dt;% 整时间节点
end