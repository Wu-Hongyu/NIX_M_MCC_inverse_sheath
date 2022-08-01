function [ mcc_para,Xsec,Xsec_H,Xsec_Hp_H2,Xsec_Hp_H,Xsec_H2p_H2,Xsec_H2p_H,Xsec_H3p_H2,Xsec_H3p_H,Xsec_Hn_H2,Xsec_Hn_H]...
                            = MCC_init( simulation ,constants,pressure,type )
%MCC_INIT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

mcc_para.T_gas=simulation.Tgas ;%K
mcc_para.pressure_gas=pressure;%Pa
mcc_para.n_target=mcc_para.pressure_gas/constants.kB/mcc_para.T_gas;
% mcc_para.n_target=1e21;
mcc_para.n_target_H=1*mcc_para.n_target;%��ԭ��1:4�ı���
mcc_para.n_target=0*mcc_para.n_target;%�����




% mcc_para.H2_Colli=1E-19;
Xsec=get_Xsec_e_H2(type, 0);
Xsec_H=get_Xsec_e_H(type, 0);

% Xsec=get_Xsec_e_H2( '2020LZS', 0);
% Xsec=get_Xsec_e_H2( 'Phelps-m', 0);
% Xsec=get_Xsec_e_H2( '1994Ness', 0);
Xsec_Hp_H2=get_Xsec_Ion('H+_H2', 0);
Xsec_Hp_H=get_Xsec_Ion('H+_H', 0);

Xsec_H2p_H2=get_Xsec_Ion('H2+_H2', 0);
Xsec_H2p_H=get_Xsec_Ion('H2+_H', 0);

% Xsec_H3p_H2=get_Xsec_Ion('H3+_H2', 0);
% Xsec_H3p_H=get_Xsec_Ion('H3+_H', 0);

% Xsec_Hn_H2=get_Xsec_Ion('H-_H2', 0);
% Xsec_Hn_H=get_Xsec_Ion('H-_H', 0);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %����Ϊ���浼��
% for i=1:length(Xsec.energy)%%����ȡ����û�����壬�޶���3000eV����
%     if Xsec.energy(i)<=3000
%         energy(i)=Xsec.energy(i);
%         tot_sec(i)=Xsec.tot{1}(i);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%���ӵ�  H2  H
[max_secXv ,~]=max(Xsec.tot{1,1}.*sqrt(2*abs(constants.q_m_ration_e)*Xsec.energy));
mcc_para.max_rate=mcc_para.n_target*max_secXv;%���Ӧ����
mcc_para.P_null=1-exp(-mcc_para.max_rate*simulation.dt);%���Ӧ����

[max_secXv_H ,~]=max(Xsec.tot{1,1}.*sqrt(2*abs(constants.q_m_ration_e)*Xsec.energy));
mcc_para.max_rate_H=mcc_para.n_target_H*max_secXv_H;%���Ӧ����
mcc_para.P_null_H=1-exp(-mcc_para.max_rate_H*simulation.dt);%���Ӧ����

%%%%%%%%%%%%%%%%%%
%%���ӵ� 
%%Hp+H2
[max_secXv_Hp_H2 ,~]=max(Xsec_Hp_H2.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_Hp)*Xsec_Hp_H2.energy));
mcc_para.max_rate_Hp_H2=mcc_para.n_target*max_secXv_Hp_H2;%���Ӧ����
mcc_para.P_null_Hp_H2=1-exp(-mcc_para.max_rate_Hp_H2*simulation.dt);%���Ӧ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Hp+H
[max_secXv_Hp_H ,~]=max(Xsec_Hp_H.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_Hp)*Xsec_Hp_H.energy));
mcc_para.max_rate_Hp_H=mcc_para.n_target_H*max_secXv_Hp_H;%���Ӧ����
mcc_para.P_null_Hp_H=1-exp(-mcc_para.max_rate_Hp_H*simulation.dt);%���Ӧ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%H2p+H2
[max_secXv_H2p_H2 ,~]=max(Xsec_H2p_H2.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_H2p)*Xsec_H2p_H2.energy));
mcc_para.max_rate_H2p_H2=mcc_para.n_target*max_secXv_H2p_H2;%���Ӧ����
mcc_para.P_null_H2p_H2=1-exp(-mcc_para.max_rate_H2p_H2*simulation.dt);%���Ӧ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%H2p+H
[max_secXv_H2p_H ,~]=max(Xsec_H2p_H.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_H2p)*Xsec_H2p_H.energy));
mcc_para.max_rate_H2p_H=mcc_para.n_target_H*max_secXv_H2p_H;%���Ӧ����
mcc_para.P_null_H2p_H=1-exp(-mcc_para.max_rate_H2p_H*simulation.dt);%���Ӧ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%H3p+H2
% [max_secXv_H3p_H2 ,~]=max(Xsec_H3p_H2.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_H3p)*Xsec_H3p_H2.energy));
% mcc_para.max_rate_H3p_H2=mcc_para.n_target*max_secXv_H3p_H2;%���Ӧ����
% mcc_para.P_null_H3p_H2=1-exp(-mcc_para.max_rate_H3p_H2*simulation.dt);%���Ӧ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%H3p+H
% [max_secXv_H3p_H ,~]=max(Xsec_H3p_H.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_H3p)*Xsec_H3p_H.energy));
% mcc_para.max_rate_H3p_H=mcc_para.n_target_H*max_secXv_H3p_H;%���Ӧ����
% mcc_para.P_null_H3p_H=1-exp(-mcc_para.max_rate_H3p_H*simulation.dt);%���Ӧ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%Hn+H2
% [max_secXv_Hn_H2 ,~]=max(Xsec_Hn_H2.tot{1,1}*sqrt(2*abs(constants.q_m_ratio_Hn)*Xsec_Hn_H2.energy));
% mcc_para.max_rate_Hn_H2=mcc_para.n_target*max_secXv_Hn_H2;%���Ӧ����
% mcc_para.P_null_Hn_H2=1-exp(-mcc_para.max_rate_Hn_H2*simulation.dt);%���Ӧ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Hp+H
% [max_secXv_Hn_H ,~]=max(Xsec_Hn_H.tot{1,1}.*sqrt(2*abs(constants.q_m_ratio_Hn)*Xsec_Hn_H.energy));
% mcc_para.max_rate_Hn_H=mcc_para.n_target_H*max_secXv_Hn_H;%���Ӧ����
% mcc_para.P_null_Hn_H=1-exp(-mcc_para.max_rate_Hn_H*simulation.dt);%���Ӧ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%e_H2
mcc_para.ion_index=length(Xsec.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index=length(Xsec.att{1,1})+mcc_para.ion_index;%���ڼ���������
mcc_para.ela_index=length(Xsec.ela{1,1})+mcc_para.att_index;%���ڼ����ǵ���
mcc_para.exc_index=length(Xsec.exc{1,1})+mcc_para.ela_index;%���ڼ����Ǽ���
%e_H
mcc_para.ion_index_H=length(Xsec_H.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index_H=length(Xsec_H.att{1,1})+mcc_para.ion_index_H;%���ڼ���������
mcc_para.ela_index_H=length(Xsec_H.ela{1,1})+mcc_para.att_index_H;%���ڼ����ǵ���
mcc_para.exc_index_H=length(Xsec_H.exc{1,1})+mcc_para.ela_index_H;%���ڼ����Ǽ���

%%%%%%%%%%%
%Hp_H2
mcc_para.ion_index_Hp_H2=length(Xsec_Hp_H2.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index_Hp_H2=length(Xsec_Hp_H2.att{1,1})+mcc_para.ion_index_Hp_H2;%���ڼ���������
mcc_para.ela_index_Hp_H2=length(Xsec_Hp_H2.ela{1,1})+mcc_para.att_index_Hp_H2;%���ڼ����ǵ���
mcc_para.exc_index_Hp_H2=length(Xsec_Hp_H2.exc{1,1})+mcc_para.ela_index_Hp_H2;%���ڼ����Ǽ���
%Hp_H
mcc_para.ion_index_Hp_H=length(Xsec_Hp_H.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index_Hp_H=length(Xsec_Hp_H.att{1,1})+mcc_para.ion_index_Hp_H;%���ڼ���������
mcc_para.ela_index_Hp_H=length(Xsec_Hp_H.ela{1,1})+mcc_para.att_index_Hp_H;%���ڼ����ǵ���
mcc_para.exc_index_Hp_H=length(Xsec_Hp_H.exc{1,1})+mcc_para.ela_index_Hp_H;%���ڼ����Ǽ���


%%%%%%%%%%%
%H2p_H2
mcc_para.ion_index_H2p_H2=length(Xsec_H2p_H2.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index_H2p_H2=length(Xsec_H2p_H2.att{1,1})+mcc_para.ion_index_H2p_H2;%���ڼ���������
mcc_para.ela_index_H2p_H2=length(Xsec_H2p_H2.ela{1,1})+mcc_para.att_index_H2p_H2;%���ڼ����ǵ���
mcc_para.exc_index_H2p_H2=length(Xsec_H2p_H2.exc{1,1})+mcc_para.ela_index_H2p_H2;%���ڼ����Ǽ���
%H2p_H
mcc_para.ion_index_H2p_H=length(Xsec_H2p_H.ion{1,1});%���ڼ����ǵ���
mcc_para.att_index_H2p_H=length(Xsec_H2p_H.att{1,1})+mcc_para.ion_index_H2p_H;%���ڼ���������
mcc_para.ela_index_H2p_H=length(Xsec_H2p_H.ela{1,1})+mcc_para.att_index_H2p_H;%���ڼ����ǵ���
mcc_para.exc_index_H2p_H=length(Xsec_H2p_H.exc{1,1})+mcc_para.ela_index_H2p_H;%���ڼ����Ǽ���

%%%%%%%%%%%
% %H3p_H2
% mcc_para.ion_index_H3p_H2=length(Xsec_H3p_H2.ion{1,1});%���ڼ����ǵ���
% mcc_para.att_index_H3p_H2=length(Xsec_H3p_H2.att{1,1})+mcc_para.ion_index_H3p_H2;%���ڼ���������
% mcc_para.ela_index_H3p_H2=length(Xsec_H3p_H2.ela{1,1})+mcc_para.att_index_H3p_H2;%���ڼ����ǵ���
% mcc_para.exc_index_H3p_H2=length(Xsec_H3p_H2.exc{1,1})+mcc_para.ela_index_H3p_H2;%���ڼ����Ǽ���
% %H3p_H
% mcc_para.ion_index_H3p_H=length(Xsec_H3p_H.ion{1,1});%���ڼ����ǵ���
% mcc_para.att_index_H3p_H=length(Xsec_H3p_H.att{1,1})+mcc_para.ion_index_H3p_H;%���ڼ���������
% mcc_para.ela_index_H3p_H=length(Xsec_H3p_H.ela{1,1})+mcc_para.att_index_H3p_H;%���ڼ����ǵ���
% mcc_para.exc_index_H3p_H=length(Xsec_H3p_H.exc{1,1})+mcc_para.ela_index_H3p_H;%���ڼ����Ǽ���

Xsec_H3p_H2=0;
Xsec_H3p_H=0;
Xsec_Hn_H2=0;
Xsec_Hn_H=0;

end

