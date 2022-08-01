function [ v_init ] = get_v_init( vth, type_init, dimension_vec )
% 初始化粒子速度
switch type_init
    case 'Maxwellian velocity'
        %按照Maxwell分布生成电子的速度分布 
        v_init=normrnd(0,vth,dimension_vec);  
    case 'Maxwellian Flux'
        error('Not Done');
        % v_init=
end
end