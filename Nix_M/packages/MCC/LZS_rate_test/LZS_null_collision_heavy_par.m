function [ collision_index , collision_type ] = LZS_null_collision_heavy_par( ve,constants,mcc_para,Xsec,vH2)
%NULL-COLLISION �˴���ʾ�йش˺�����ժҪ
%���ط��������ĵ��ӵ�����

% select_sec_tot=interp1(energy,tot_sec,select_energy,'Linear');
% select_sec_ela=interp1(Xsec.ela{1,1}{1,1}(:,1),Xsec.ela{1,1}{1,1}(:,2),select_Te,'Linear');���󣬲�������ѡ��������ʵĵ����¶ȵ�ֵ���м��㣬�õ������ѡȡ�ĵ��ӵ��������㡣

% % P_ela=select_sec_ela/select_sec_tot;���󣬲������ý���ı�ֵ���������ʵı�ֵ
max_collision_num=round(mcc_para.P_null*length(ve));
%��������ײ��

null_index=randperm(length(ve),max_collision_num);%��ȡ�����˿���ײ������

null_H2_index=randperm(length(vH2),max_collision_num);%��ȡ�����˿���ײ������

v_targetFrame=ve(null_index,:)-vH2(null_H2_index,:);

EineV=@(v,q_m_ratio) (v(:,1).^2+v(:,2).^2+v(:,3).^2)/(2*abs(q_m_ratio));%�ٶȱ�����

select_energy=EineV(v_targetFrame,2*constants.q_m_ration_H);%��ȡ�����ӵ�����ֵ

% select_energy=50*ones(length(select_energy),1);

for i=1:length(Xsec.ion{1,1})
    select_sec_ion(:,i)=interp1(Xsec.ion{1,1}{1,i}(:,1),Xsec.ion{1,1}{1,i}(:,2),select_energy,'Linear');%��ȡ�����ӵ�����ֵ��Ӧ�Ľ����ֵ
    rate_ion(:,i)=mcc_para.n_target*select_sec_ion(:,i).*sqrt(2*abs(constants.q_m_ration_H*2)*select_energy);%��ȡ�����ӵĵ���Ƶ��
end

for i=1:length(Xsec.att{1,1})
    select_sec_att(:,i)=interp1(Xsec.att{1,1}{1,i}(:,1),Xsec.att{1,1}{1,i}(:,2),select_energy,'Linear');%��ȡ�����ӵ�����ֵ��Ӧ�Ľ����ֵ
    rate_att(:,i)=mcc_para.n_target*select_sec_att(:,i).*sqrt(2*abs(constants.q_m_ration_H*2)*select_energy);%��ȡ�����ӵĵ���Ƶ��
end

for i=1:length(Xsec.ela{1,1})
    select_sec_ela(:,i)=interp1(Xsec.ela{1,1}{1,i}(:,1),Xsec.ela{1,1}{1,i}(:,2),select_energy,'Linear');%��ȡ�����ӵ�����ֵ��Ӧ�Ľ����ֵ
    rate_ela(:,i)=mcc_para.n_target*select_sec_ela(:,i).*sqrt(2*abs(constants.q_m_ration_H*2)*select_energy);%��ȡ�����ӵĵ���Ƶ��
end

% yy=arrayfun(@(i)interp1(Xsec.exc{1,1}{1,i}(:,1),Xsec.exc{1,1}{1,i}(:,2),select_energy(1,1)),1:length(Xsec.exc{1,1}));
% for i=1:length(Xsec.exc{1,1})
% cross_sec(:,i)=Xsec.exc{1,1}{1,i}(:,2);
% sec_energy(:,i)=Xsec.exc{1,1}{1,i}(:,1);
% end

for i=1:length(Xsec.exc{1,1})
    select_sec_exc(:,i)=interp1(Xsec.exc{1,1}{1,i}(:,1),Xsec.exc{1,1}{1,i}(:,2),select_energy,'Linear');%��ȡ�����ӵ�����ֵ��Ӧ�Ľ����ֵ
    rate_exc(:,i)=mcc_para.n_target*select_sec_exc(:,i).*sqrt(2*abs(constants.q_m_ration_H*2)*select_energy);%��ȡ�����ӵĵ���Ƶ��
end
%%%%%%%%%%%������ʵ��������ֵ�ж�һ���Ƿ��ܷ�Ӧ�����ܷ�Ӧ�ľͲ����в�ֵ�Ĵ�����

% select_sec_ela2=interp1(Xsec.ela{1,1}{1,1}(:,1),Xsec.ela{1,1}{1,1}(:,2),select_energy,'Linear');%��ȡ�����ӵ�����ֵ��Ӧ�Ľ����ֵ
% rate_ela2=mcc_para.n_target*select_sec_ela.*sqrt(2*abs(constants.q_m_ration_e)*select_energy);%��ȡ�����ӵĵ���Ƶ��
P_ion=rate_ion./mcc_para.max_rate;%�������
P_att=rate_att./mcc_para.max_rate;
P_ela=rate_ela./mcc_para.max_rate;
P_exc=rate_exc./mcc_para.max_rate;

% P_ion=0.1*ones(length(P_ion),1);
% P_att=0.09*ones(length(P_att),1);
% P_ela=0.08*ones(length(P_ela),1);
% for jj=1:15
% P_exc(:,jj)=0.07*((16-jj)/15)*ones(length(P_exc(:,jj)),1);
% end





P_accum=cumsum([P_ion P_att P_ela P_exc],2); %�����ۼ�
% ela_index=null_index(rand(max_collision_num,1)<P_ela);%ȡ�����С��v_ela/v_tot�����ӵ�����
random_number=rand(max_collision_num,1);

reaction_type=[];
for n=1:max_collision_num
    for j=1:length(P_accum(1,:))
        if random_number(n)<P_accum(n,j)%�ж��Ƿ�����Ӧ���Լ�������ʲô��Ӧ
            reaction_type(n,1)=j;
            break
        end
    end
end

if isempty(reaction_type)
    collision_index=[];
    collision_type=[];
else
    happen=logical(reaction_type);%��ѡ�����˷�Ӧ������
    collision_index=null_index(happen);%��Ӧ��������ve�е�index
    collision_type=reaction_type(happen);%��Ӧ�ķ�Ӧ����  1ion    2att   3ela   4exc
end



