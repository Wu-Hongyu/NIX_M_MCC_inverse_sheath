function [accur_post_v] = Ion_ela_process(vp)
%ION_ELA_PROCESS 此处显示有关此函数的摘要
%   此处显示详细说明
if ~isempty(vp)
    v_value=sqrt(vp(:,1).^2+vp(:,2).^2+vp(:,3).^2);%速度模值
    %     EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
    %     switch react
    %         case 'Hp'
    %             E=EineV(vp(:,:),constants.q_m_ratio_Hp);%发生弹性碰撞的粒子的入射能量
    %         case 'H2p'
    %             E=EineV(vp(:,:),constants.q_m_ratio_H2p);%发生弹性碰撞的粒子的入射能量
    %     end
    alpha=rand(length(vp(:,1)),1);
    cost=sqrt(alpha);
    sint=sqrt(1-cost.^2);
    p=2*pi*rand(length(vp(:,1)),1);%φ角随机0~2pi
    cosp=cos(p);
    sinp=sin(p);
    
    v_vertical=sqrt(vp(:,1).^2+vp(:,2).^2);%王辉辉2-36  2-37
    dvx=vp(:,1)./v_vertical.*vp(:,3).*sint.*cosp-vp(:,2)./v_vertical.*v_value.*sint.*sinp-vp(:,1).*(1-cost);
    dvy=(vp(:,2)./v_vertical.*vp(:,3).*sint.*cosp+vp(:,1)./v_vertical.*v_value.*sint.*sinp-vp(:,2).*(1-cost));
    dvz=-v_vertical.*sint.*cosp-vp(:,3).*(1-cost);
    accur_post_v=(vp(:,:)+[dvx dvy dvz]).*cost;
else
    accur_post_v=[];
end
end

