function [ v_init ] = get_v_init( vth, type_init, dimension_vec )
% ��ʼ�������ٶ�
switch type_init
    case 'Maxwellian velocity'
        %����Maxwell�ֲ����ɵ��ӵ��ٶȷֲ� 
        v_init=normrnd(0,vth,dimension_vec);  
    case 'Maxwellian Flux'
        error('Not Done');
        % v_init=
end
end