function [ accur_post_v ] = exc_process( exc_index,exc_react_order,ve,constants ,  excThresh )
%EXC_PROCESS 此处显示有关此函数的摘要
%   此处显示详细说明
accur_post_v=[];
if ~isempty(exc_index)

% z=[zeros(size(exc_index,2),2) ones(size(exc_index,2),1) ];
EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量
E_loss=[excThresh{1,1}{1,exc_react_order}]';
E=EineV(ve(exc_index,:),constants.q_m_ratio_e);

E_remain=E-E_loss;

v_inc=ve(exc_index,:).*sqrt(E_remain./E);
cost=(2+E_remain-2*(1+E_remain).^rand(length(exc_index),1))./E_remain;
sint=sqrt(1-cost.^2);
v_value=sqrt(v_inc(:,1).^2+v_inc(:,2).^2+v_inc(:,3).^2);

phi=2*pi*rand(length(exc_index),1);
cosp=cos(phi);
sinp=sin(phi);
v_vertical=sqrt(v_inc(:,1).^2+v_inc(:,2).^2);
dvx=v_inc(:,1)./v_vertical.*v_inc(:,3).*sint.*cosp-v_inc(:,2)./v_vertical.*v_value.*sint.*sinp-v_inc(:,1).*(1-cost);
dvy=(v_inc(:,2)./v_vertical.*v_inc(:,3).*sint.*cosp+v_inc(:,1)./v_vertical.*v_value.*sint.*sinp-v_inc(:,2).*(1-cost));
dvz=-v_vertical.*sint.*cosp-v_inc(:,3).*(1-cost);
% dEtoE=constants.me/(constants.mH+constants.me)*(1-cost);
accur_post_v=(v_inc(:,:)+[dvx dvy dvz]);

% E3=EineV(accur_post_v,constants.q_m_ration_e);

end


end

