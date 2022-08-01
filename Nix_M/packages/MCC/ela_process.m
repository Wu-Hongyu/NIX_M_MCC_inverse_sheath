function [ accur_post_v ] = ela_process( ela_index,ve,constants,react)
%ELA_ 此处显示有关此函数的摘要
%   此处显示详细说明
if ~isempty(ela_index)
    
    v_value=sqrt(ve(ela_index,1).^2+ve(ela_index,2).^2+ve(ela_index,3).^2);%速度模值
    
    % z=[zeros(size(ela_index,2),2) ones(size(ela_index,2),1) ];
    
    EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
    E=EineV(ve(ela_index,:),constants.q_m_ratio_e);%发生弹性碰撞的粒子的入射能量
    
    cost=(2+E-2*(1+E).^rand(length(ela_index),1))./E;%散射角采样%V.  Vahedi,  M.  Surendra.  A  Monte  Carlo  collision  model  for  the  particle-in-cell  method: applications to argon and oxygen discharges
    %     cost=cos(pi*rand(length(ela_index),1));
    sint=sqrt(1-cost.^2);
    
    % p=2*pi*rand(length(ela_index),1);
    % cosp=cos(p);
    % sinp=sin(p);
    % v_direction=ve(ela_index,:)./v_value;
    % costz=v_direction(:,3);
    % sintz=sqrt(1-costz.^2);
    % rough_post_direction=v_value.*(cross(v_direction,z,2).*sint.*sinp./sintz+cross(v_direction,cross(z,v_direction,2),2).*sint.*cosp./sintz+v_direction.*cost);
    %     E2=EineV(rough_post_direction,constants.q_m_ration_e);%暂时不考虑弹碰传递给中性粒子的能量损失，
    
    % p=p+pi;%上下两种方法定义的角度phi相差了180°
    p=2*pi*rand(length(ela_index),1);%φ角随机0~2pi
    cosp=cos(p);
    sinp=sin(p);
    v_vertical=sqrt(ve(ela_index,1).^2+ve(ela_index,2).^2);%王辉辉2-36  2-37
    dvx=ve(ela_index,1)./v_vertical.*ve(ela_index,3).*sint.*cosp-ve(ela_index,2)./v_vertical.*v_value.*sint.*sinp-ve(ela_index,1).*(1-cost);
    dvy=(ve(ela_index,2)./v_vertical.*ve(ela_index,3).*sint.*cosp+ve(ela_index,1)./v_vertical.*v_value.*sint.*sinp-ve(ela_index,2).*(1-cost));
    dvz=-v_vertical.*sint.*cosp-ve(ela_index,3).*(1-cost);
    switch react
        case  'H2'
            dEtoE=constants.me/(constants.mH2+constants.me)*(1-cost);%王辉辉3-2 弹性碰撞中的电子能量损失值比例
        case 'H'
            dEtoE=constants.me/(constants.mH+constants.me)*(1-cost);%王辉辉3-2 弹性碰撞中的电子能量损失值比例
    end
%     accur_post_v=(ve(ela_index,:)+[dvx dvy dvz]).*(1-dEtoE);
    accur_post_v=(ve(ela_index,:)+[dvx dvy dvz]);
    %   E3=EineV(accur_post_v,constants.q_m_ration_e);
    
    
else
    accur_post_v=[];
end

end
