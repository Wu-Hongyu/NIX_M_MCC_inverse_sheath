function [ ion_inc_post_v ,ion_gen_post_v ] = ion_process( ion_index,ion_react_order ,ve,constants ,  ionThresh )
%ION_PROCESS 此处显示有关此函数的摘要
%   此处显示详细说明

ion_inc_post_v=[];
ion_gen_post_v=[];
if ~isempty(ion_index)
    
    % z=[zeros(size(ion_index,2),2) ones(size(ion_index,2),1) ];
    EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
    E=EineV(ve(ion_index,:),constants.q_m_ratio_e);
    E_loss=[ionThresh{1,1}{1,ion_react_order}]';%不同的反应阈值
    E_remain=E-E_loss;
    
    E_inc=(E_remain).*rand(length(ion_index),1);
    E_gen=E_remain-E_inc;%剩余能量随机分配给入射和电离电子
    E_all=[E_inc; E_gen];
    
    v_inc=ve(ion_index,:).*sqrt(E_inc./E);
    v_gen=ve(ion_index,:).*sqrt(E_gen./E);
    v_all=[v_inc; v_gen];
    cost=(2+E_all-2*(1+E_all).^rand(length(E_all),1))./E_all;
    sint=sqrt(1-cost.^2);
    v_value=sqrt(v_all(:,1).^2+v_all(:,2).^2+v_all(:,3).^2);
    
    phi=2*pi*rand(length(E_all),1);
    cosp=cos(phi);
    sinp=sin(phi);
    v_vertical=sqrt(v_all(:,1).^2+v_all(:,2).^2);
    dvx=v_all(:,1)./v_vertical.*v_all(:,3).*sint.*cosp-v_all(:,2)./v_vertical.*v_value.*sint.*sinp-v_all(:,1).*(1-cost);
    dvy=(v_all(:,2)./v_vertical.*v_all(:,3).*sint.*cosp+v_all(:,1)./v_vertical.*v_value.*sint.*sinp-v_all(:,2).*(1-cost));
    dvz=-v_vertical.*sint.*cosp-v_all(:,3).*(1-cost);
%     dEtoE=constants.me/(constants.mH+constants.me)*(1-cost);
    ion_all_post_v=(v_all(:,:)+[dvx dvy dvz]);
    ion_inc_post_v=ion_all_post_v(1:length(E_inc),:);
    ion_gen_post_v=ion_all_post_v(length(E_inc)+1:end,:);
    % E3=EineV(ion_all_post_v,constants.q_m_ration_e);
    
    
end


end

