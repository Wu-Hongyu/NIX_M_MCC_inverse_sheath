function [ accur_post_v ] = ela_process( ela_index,ve,constants,react)
%ELA_ �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if ~isempty(ela_index)
    
    v_value=sqrt(ve(ela_index,1).^2+ve(ela_index,2).^2+ve(ela_index,3).^2);%�ٶ�ģֵ
    
    % z=[zeros(size(ela_index,2),2) ones(size(ela_index,2),1) ];
    
    EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%�ٶȱ�����
    E=EineV(ve(ela_index,:),constants.q_m_ratio_e);%����������ײ�����ӵ���������
    
    cost=(2+E-2*(1+E).^rand(length(ela_index),1))./E;%ɢ��ǲ���%V.  Vahedi,  M.  Surendra.  A  Monte  Carlo  collision  model  for  the  particle-in-cell  method: applications to argon and oxygen discharges
    %     cost=cos(pi*rand(length(ela_index),1));
    sint=sqrt(1-cost.^2);
    
    % p=2*pi*rand(length(ela_index),1);
    % cosp=cos(p);
    % sinp=sin(p);
    % v_direction=ve(ela_index,:)./v_value;
    % costz=v_direction(:,3);
    % sintz=sqrt(1-costz.^2);
    % rough_post_direction=v_value.*(cross(v_direction,z,2).*sint.*sinp./sintz+cross(v_direction,cross(z,v_direction,2),2).*sint.*cosp./sintz+v_direction.*cost);
    %     E2=EineV(rough_post_direction,constants.q_m_ration_e);%��ʱ�����ǵ������ݸ��������ӵ�������ʧ��
    
    % p=p+pi;%�������ַ�������ĽǶ�phi�����180��
    p=2*pi*rand(length(ela_index),1);%�ս����0~2pi
    cosp=cos(p);
    sinp=sin(p);
    v_vertical=sqrt(ve(ela_index,1).^2+ve(ela_index,2).^2);%���Ի�2-36  2-37
    dvx=ve(ela_index,1)./v_vertical.*ve(ela_index,3).*sint.*cosp-ve(ela_index,2)./v_vertical.*v_value.*sint.*sinp-ve(ela_index,1).*(1-cost);
    dvy=(ve(ela_index,2)./v_vertical.*ve(ela_index,3).*sint.*cosp+ve(ela_index,1)./v_vertical.*v_value.*sint.*sinp-ve(ela_index,2).*(1-cost));
    dvz=-v_vertical.*sint.*cosp-ve(ela_index,3).*(1-cost);
    switch react
        case  'H2'
            dEtoE=constants.me/(constants.mH2+constants.me)*(1-cost);%���Ի�3-2 ������ײ�еĵ���������ʧֵ����
        case 'H'
            dEtoE=constants.me/(constants.mH+constants.me)*(1-cost);%���Ի�3-2 ������ײ�еĵ���������ʧֵ����
    end
%     accur_post_v=(ve(ela_index,:)+[dvx dvy dvz]).*(1-dEtoE);
    accur_post_v=(ve(ela_index,:)+[dvx dvy dvz]);
    %   E3=EineV(accur_post_v,constants.q_m_ration_e);
    
    
else
    accur_post_v=[];
end

end
