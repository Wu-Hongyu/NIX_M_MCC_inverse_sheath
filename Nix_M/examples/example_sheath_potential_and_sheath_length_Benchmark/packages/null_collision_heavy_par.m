function [ collision_index, collision_type ] = null_collision_heavy_par( vP,constants,mcc_para,Xsec,react)
%NULL-COLLISION 此处显示有关此函数的摘要
%返回发生弹碰的电子的索引
switch react
    % select_sec_tot=interp1(energy,tot_sec,select_energy,'Linear');
    % select_sec_ela=interp1(Xsec.ela{1,1}{1,1}(:,1),Xsec.ela{1,1}{1,1}(:,2),select_Te,'Linear');错误，并不是用选定最大速率的电子温度的值进行计算，用的是随机选取的电子的能量计算。
    case 'Hp+H2'
        P_null=mcc_para.P_null_Hp_H2;
        q_m_ratio=constants.q_m_ratio_Hp;
        target_n=mcc_para.n_target;
        max_rate=mcc_para.max_rate_Hp_H2;
    case 'Hp+H'
        P_null=mcc_para.P_null_Hp_H;
        q_m_ratio=constants.q_m_ratio_Hp;
        target_n=mcc_para.n_target_H;
        max_rate=mcc_para.max_rate_Hp_H;
    case 'H2p+H2'
        P_null=mcc_para.P_null_H2p_H2;
        q_m_ratio=constants.q_m_ratio_H2p;
        target_n=mcc_para.n_target;
        max_rate=mcc_para.max_rate_H2p_H2;
    case 'H2p+H'
        P_null=mcc_para.P_null_H2p_H;
        q_m_ratio=constants.q_m_ratio_H2p;
        target_n=mcc_para.n_target_H;
        max_rate=mcc_para.max_rate_H2p_H;
    case 'H3p+H2'
        P_null=mcc_para.P_null_H3p_H2;
        q_m_ratio=constants.q_m_ratio_H3p;
        target_n=mcc_para.n_target;
        max_rate=mcc_para.max_rate_H3p_H2;
    case 'H3p+H'
        P_null=mcc_para.P_null_H3p_H;
        q_m_ratio=constants.q_m_ratio_H3p;
        target_n=mcc_para.n_target_H;
        max_rate=mcc_para.max_rate_H3p_H;
    case 'Hn+H2'
        P_null=mcc_para.P_null_Hn_H2;
        q_m_ratio=constants.q_m_ratio_Hp;
        target_n=mcc_para.n_target;
        max_rate=mcc_para.max_rate_Hn_H2;
    case 'Hn+H'
        P_null=mcc_para.P_null_Hn_H;
        q_m_ratio=constants.q_m_ratio_Hp;
        target_n=mcc_para.n_target_H;
        max_rate=mcc_para.max_rate_Hn_H;
end
% % P_ela=select_sec_ela/select_sec_tot;错误，并不是用截面的比值，是用速率的比值
max_collision_num=round(P_null*length(vP));
%最大可能碰撞数

null_index=randperm(length(vP),max_collision_num);%抽取发生了空碰撞的粒子

% null_H2_index=randperm(length(vH2),max_collision_num);%抽取发生了空碰撞的粒子

v_targetFrame=vP(null_index,:);

EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%速度变能量

select_energy=EineV(v_targetFrame,q_m_ratio);%抽取的粒子的能量值

% select_energy=50*ones(length(select_energy),1);

% % % %         rate_colli=mcc_para.n_target*mcc_para.H2_Colli*sqrt(2*abs(q_m_ratio)*select_energy);%抽取的粒子的弹碰频率


for i=1:length(Xsec.ion{1,1})
    select_sec_ion(:,i)=interp1(Xsec.ion{1,1}{1,i}(:,1),Xsec.ion{1,1}{1,i}(:,2),select_energy,'Linear');%抽取的粒子的能量值对应的截面插值
    rate_ion(:,i)=target_n*select_sec_ion(:,i).*sqrt(2*abs(q_m_ratio)*select_energy);%抽取的粒子的弹碰频率
end

for i=1:length(Xsec.att{1,1})
    select_sec_att(:,i)=interp1(Xsec.att{1,1}{1,i}(:,1),Xsec.att{1,1}{1,i}(:,2),select_energy,'Linear');%抽取的粒子的能量值对应的截面插值
    rate_att(:,i)=target_n*select_sec_att(:,i).*sqrt(2*abs(q_m_ratio)*select_energy);%抽取的粒子的弹碰频率
end

for i=1:length(Xsec.ela{1,1})
    select_sec_ela(:,i)=interp1(Xsec.ela{1,1}{1,i}(:,1),Xsec.ela{1,1}{1,i}(:,2),select_energy,'Linear');%抽取的粒子的能量值对应的截面插值
    rate_ela(:,i)=target_n*select_sec_ela(:,i).*sqrt(2*abs(q_m_ratio)*select_energy);%抽取的粒子的弹碰频率
end

for i=1:length(Xsec.exc{1,1})
    select_sec_exc(:,i)=interp1(Xsec.exc{1,1}{1,i}(:,1),Xsec.exc{1,1}{1,i}(:,2),select_energy,'Linear');%抽取的粒子的能量值对应的截面插值
    rate_exc(:,i)=target_n*select_sec_exc(:,i).*sqrt(2*abs(q_m_ratio)*select_energy);%抽取的粒子的弹碰频率
end

%%%%%%%%%%%这里其实可以用阈值判断一下是否能反应，不能反应的就不进行插值的处理了

P_ion=rate_ion./max_rate;%电离概率
P_att=rate_att./max_rate;
P_ela=rate_ela./max_rate;
P_exc=rate_exc./max_rate;

P_accum=cumsum([P_ion P_att P_ela P_exc],2); %概率累加
% select_sec_ela2=interp1(Xsec.ela{1,1}{1,1}(:,1),Xsec.ela{1,1}{1,1}(:,2),select_energy,'Linear');%抽取的粒子的能量值对应的截面插值
% rate_ela2=mcc_para.n_target*select_sec_ela.*sqrt(2*abs(constants.q_m_ration_e)*select_energy);%抽取的粒子的弹碰频率
%         P_Colli=rate_colli./mcc_para.max_rate_H2Colli;%电离概率

% ela_index=null_index(rand(max_collision_num,1)<P_ela);%取随机数小于v_ela/v_tot的粒子的索引
random_number=rand(max_collision_num,1);

reaction_type=[];
for n=1:max_collision_num
    for j=1:length(P_accum(1,:))
        if random_number(n)<P_accum(n,j)%判断是否发生反应，以及发生了什么反应
            reaction_type(n,1)=j;
            break
        end
    end
end

if isempty(reaction_type)
    collision_index=[];
    collision_type=[];
else
    happen=logical(reaction_type);%挑选发生了反应的粒子
    collision_index=null_index(happen);%反应的粒子在ve中的index
    collision_type=reaction_type(happen);%对应的反应类型  1ion    2att   3ela   4exc
end

%         colli_happen=(random_number<P_Colli);
%         collision_index=[];
%         if sum(colli_happen)~=0
%             collision_index=null_index(colli_happen);
%         end

% end