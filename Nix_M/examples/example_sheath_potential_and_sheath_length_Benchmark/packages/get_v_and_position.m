function [ v, position ] = get_v_and_position( v,position,q_m_ratio, type,E_particle,B_particle, simulation )
% particle pusher的package，各种pusher函数的集中接口
% 传值性能低。仅在不在意性能时使用，如测试功能时
switch type
    case '1D leap-frog' %1D仅电场
        [ v, position ] = pusher1D_leap_frog( v, position, q_m_ratio, E_particle ,simulation );
    case '1D leap-frog+B by vector analysis'
        [ v, position ] = pusher1D_leap_frog_B_vec( v, position, q_m_ratio, E_particle,B_particle,simulation );
    case 'parallel'
        error('Not Done');
    case 'Boris' % 2D
        error('Not Done');
end
end