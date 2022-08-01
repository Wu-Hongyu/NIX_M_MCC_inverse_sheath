function [ Xsec ] = get_Xsec_e_H2( type, flag_display)
% ����ImportLXCat�õ���e-H2ԭʼ���ݣ���ȡMCC�������
% ʹ��
% type��e_H2ɢ�����������������һ�㼴���ݿ���
% flag_display��0�򲻻�ͼ�������Ӧ���࣬1����ʾ
% ע�⣬���ļ�����.\others\choose_Xsec_e_n_of_hydrogen_plasma
% \compare_Xsec_e_H2.mʵ��
% �����ԵıȽ���compare_Xsec_e_H2
% �����Ե�test��test_get_Xsec

% ͨ��ȫ���滻E_head��E_endֵ�����޸Ľ����������
path_list=get_path_init('test000');
recation_type={['e-H2-' type]};  %�����cell
switch type       
    case '2020LZS' % ����ɽ
        % ����ʹ��CCC
        % ������ת����
        % �񶯼���ʹ��2003Janev��ѡȡv1,v2
        % ���Ӽ���ʹ��CCC�����ݹ��̵�N�ͽ���ѡ��13�����̣�
        % LZS����Fantz2006������ֵ
        % ����ʹ��2011Wuenderlich��ѡȡ�ǽ����Ե��롣����չ��������
        % ����ʹ��2003Janev��ѡȡv=0������չ��v=0~14������Eth�������E_head
        
        %%%%%%% ����CCCԭʼ����
        data_dir={[path_list.cross_sections '\e_H2_CCC_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% ׼������ԭʼ����
        idx_use=[]; % ���ù������
        E_head=1e-3; % ����������Χ
        E_end=1e4; 
        
        %%%%%%% ���Դ���
        for E_end_i=logspace(log10(400),log10(E_end),10)
            Xsec=Xsec.add_key_points2(E_head, E_end_i);
        end
        % CCC����ֻ��300eV����ֻ��E_end����һ���㣬trapzʱ���ܺ��Ե���
        
        %%%%%%% �񶯼�������
        % ����excĩβ����ѡ��ʱ��ǰ
        num0_exc=length(Xsec.exc{1});
        % ������ϱ��ʽ 2003Janev-ch 4.1.1-eq.74
        xsec_v=@(ksv,dE,E) 5.78*ksv/(dE^4)*((dE/E)^2)...
            *(1-dE/E)^(6.11/sqrt(dE))*1E-16*1E-4;
        % dE: ��E��the energy difference of v = 0 and v'
        % ksv: factor
        % v=0->v=1
        Eth_v1=0.516; % �Ԧ�EΪ��ֵ
        xsec_v1=@(E) xsec_v(1,Eth_v1,E);
        x_v1=logspace(log10(Eth_v1),log10(E_end),100); 
        y_v1=arrayfun(xsec_v1,x_v1);
        % ����Xsec
        k=num0_exc+1;
        Xsec.exc{1}{k}=[x_v1',y_v1'];
        Xsec.excThresh{1}{k}=Eth_v1;
        Xsec.info.exc{1}{k}='E+H2->E+H2(v=1),Excitation';
        idx_use=[idx_use k];
        % v=0->v=2
        Eth_v2=1.003; % �Ԧ�EΪ��ֵ
        xsec_v2=@(E) xsec_v(0.628,Eth_v2,E);
        x_v2=logspace(log10(Eth_v2),log10(E_end),100); 
        y_v2=arrayfun(xsec_v2,x_v2);
        % ����Xsec
        k=num0_exc+2;
        Xsec.exc{1}{k}=[x_v2',y_v2'];
        Xsec.excThresh{1}{k}=Eth_v2;
        Xsec.info.exc{1}{k}='E+H2->E+H2(v=2),Excitation';
        idx_use=[idx_use k];
        
        %%%%%%% ���Ӽ�������
        % CCC exc 1ΪGTSC���ܽ���
        % ��һ��ele��LZS���ݽ�����k��Сѡȡ��ΪCCC����ele
        % LZS��2006Fantzb����̬������Ϊ��ֵ
        Xsec.excThresh{1}{2}=11.90; %a3/11.90
        Xsec.excThresh{1}{3}=14.62; %B"1/14.62
        Xsec.excThresh{1}{4}=13.84; %B'1/13.84
        Xsec.excThresh{1}{5}=11.37; %B1/11.37
        Xsec.excThresh{1}{7}=12.41; %C1/12.41
        Xsec.excThresh{1}{8}=11.89; %c3/11.89
        Xsec.excThresh{1}{10}=14.74; %D'1/14.74
        Xsec.excThresh{1}{11}=14.13; %D1/14.13
        Xsec.excThresh{1}{12}=13.98; %d3/13.98
        Xsec.excThresh{1}{13}=13.36; %e3/13.36
        Xsec.excThresh{1}{14}=12.42; %EF1/12.42
        Xsec.excThresh{1}{17}=13.98; %h3/13.98
        idx_use=[idx_use 2:5,7:8,10:14,17];
        % �ֽ�
        % ����ײ�����µ�����ӷֽⷴӦ��һ�㾭�ɼ�����
        % ���ȷ���������Ӧ��Ȼ���ټ���̬����ӷֽ�
        % ����MCC��ģ��H2/H����˿��Բ����ǵڶ���
        % �˴����俼��b3���Ӽ�����1987Janev����ֵ
        Xsec.excThresh{1}{20}=8.5; %b3
        idx_use=[idx_use 20];

        % �ڶ���ele��PengChen����k*��Eѡȡ����CCCȫ��ele
%         % ����2003Janev������ֵ������h3��J1ʹ��2006Fantzb����̬����
%         p_array=[2:5 7:8 10:20];
%         Eth_array=[11.72, 15.47, 14.85, 12.754, 13.29, 11.72, 15.555, 14.996, 13.6, 13, 13.13, 14.816, 14.98, 13.98, 14.824, 14.07, 7.93];
%         for j=1:length(p_array)
%             Xsec.excThresh{1}{p_array(j)}=Eth_array(j);
%         end
%         idx_use=[idx_use p_array];

        % ѡ��Ӧ
        Xsec=Xsec.choose_exc(1, idx_use);
        
        for E_end_i=logspace(log10(300),log10(E_end),10)
            Xsec=Xsec.add_key_points_exc(E_head, E_end_i);
        end
        
        %%%%%%% ���봦��
        % ��LZS���ο�2011Wuenderlich
        % Ӧ����2011Wuenderlich���ۼ�������LZS��Emin��1e3eV�����
        % ����
        a10=[4.67496911419161E-24 7.48967037030240E-26 2.27812062678045E-24 ];
        a9=[-2.32941049054735E-22 -3.77513517410756E-24 -1.34333647539531E-22];
        a8=[5.18506369967717E-21 8.50625788410889E-23 3.54455852390902E-21];
        a7=[-6.79075282615963E-20 -1.12844253900045E-21 -5.51074353507411E-20];
        a6=[5.79670183018440E-19 9.76260841920836E-21 5.58984626136045E-19];
        a5=[-3.37142238452710E-18 -5.75725333009331E-20 -3.86510112579526E-18];
        a4=[1.35382511176976E-17 2.34479544100933E-19 1.84475186160132E-17];
        a3=[-3.70855112075379E-17 -6.51526822247295E-19 -6.00035106992575E-17];
        a2=[6.63506283888535E-17 1.18236733917101E-18 1.27269671523665E-16];
        a1=[-6.99974257106332E-17 -1.26530045500246E-18 -1.58914992792416E-16];
        a0=[3.30287142988073E-17 6.05859716663196E-19 8.86780631708612E-17];
        coeff_ion=[a10;a9;a8;a7;a6;a5;a4;a3;a2;a1;a0];
        Emin_ion=[16.01 18.47 41.98]; % ��Ӧ�ڷ�0���� ��Emin��
        Emax_ion=1e3;
        xsec_ion=@(E,i) sum( coeff_ion(:,i).*(log(E)).^(10:-1:0)' );
        % ����-�ǽ�����
        Eth_ion1=15.4; % LZS�ĵ��ϵ���ֵ
        xsec_ion1=@(E) xsec_ion(E,1);
        x_ion1=logspace(log10(Emin_ion(1)),log10(Emax_ion),100); 
        y_ion1=arrayfun(xsec_ion1,x_ion1);
        % ����Xsec
        k=1;
        Xsec.ion{1}{k}=[x_ion1',y_ion1'];
        Xsec.ionThresh{1}{k}=Eth_ion1;
        Xsec.info.ion{1}{k}='E+H2->2E+H2+,Ionization';
        % ����-������ g ����С��1e22������ֵ�ϸߣ�����
        Eth_ion2=18.15; % LZS�ĵ��ϵ���ֵ
        xsec_ion2=@(E) xsec_ion(E,2);
        x_ion2=logspace(log10(Emin_ion(2)),log10(Emax_ion),100); 
        y_ion2=arrayfun(xsec_ion2,x_ion2);
        % ����Xsec
        k=k+1;
        Xsec.ion{1}{k}=[x_ion2',y_ion2'];
        Xsec.ionThresh{1}{k}=Eth_ion2;
        Xsec.info.ion{1}{k}='E+H2->2e+H2+(2��g+)->2E+H+H+,Ionization';
        % ����-������ u
        Eth_ion3=30.6; % LZS�ĵ��ϵ���ֵ
        xsec_ion3=@(E) xsec_ion(E,3);
        x_ion3=logspace(log10(Emin_ion(3)),log10(Emax_ion),100); 
        y_ion3=arrayfun(xsec_ion3,x_ion3);
        % ����Xsec
        k=k+1;
        Xsec.ion{1}{k}=[x_ion3',y_ion3'];
        Xsec.ionThresh{1}{k}=Eth_ion3;
        Xsec.info.ion{1}{k}='E+H2->2e+H2+(2��u+)->2E+H+H+,Ionization';
        % ͨ������������Դ�Ƚϣ�test ok
        for E_end_i=logspace(log10(Emax_ion),log10(E_end),10)
            Xsec=Xsec.add_key_points_ion(E_head, E_end_i);
        end
        
        %%%%%%% ��������
        % ϵ�� 2003Janev-Tab 27
        Eth_att=[3.72, 3.21, 2.72, 2.26, 1.83, 1.43, 1.36, 0.713, 0.397, 0.113, -0.139, -0.354, -0.529, -0.659, -0.736];
        sigma_max_att=[0.0000322, 0.000518, 0.00416, 0.022, 0.122, 0.453, 1.51, 4.48, 10.1, 13.9, 11.8, 8.87, 7.11, 5, 3.35];
        E0_att=0.45;
        % ������ϱ��ʽ 2003Janev-ch 4.4.1-eq.124
        xsec_att=@(E,v) sigma_max_att(v+1)*1E-20*exp(-(E-abs(Eth_att(v+1)))/E0_att);
%         for v=0:14
        for v=0:0 %��ʱֻ����v=0��̬
            xsec_att_v=@(E) xsec_att(E,v);
            Eth_att_v=Eth_att(v+1);
            if Eth_att_v>0
                x_att_v=logspace(log10(Eth_att_v+1e-3),log10(E_end),100);
            else
                Eth_att_v=E_head;
                x_att_v=logspace(log10(E_head),log10(E_end),100);
            end
            y_att_v=arrayfun(xsec_att_v,x_att_v);
            % ����Xsec
            k=v+1;
            Xsec.att{1}{k}=[x_att_v',y_att_v'];
            Xsec.attThresh{1}{k}=Eth_att_v;
            Xsec.info.att{1}{k}=['E+H2(v=' num2str(v) ')->H+H-,Attachment'];
        end
        % test�����Nix_M���������׶Ա�.pptx
        
        %%%%%%% ��������
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %��info�������ֵ��Ϣ
    case 'Phelps-m'
        % Phelps�����̵�Refϸ�ڼ�pengchen�ʼ�
        % ����ʹ��IST��Phelps��ela��eff���㣬��bug��
        % ��ת����ʹ��Phelps����j0-2��j1-3
        % �񶯼���ʹ��Phelps��ѡ��v1,2,3
        % ���Ӽ���ʹ��Phelps��ѡ��b3,B1,c3,a3,C1,D3,N��4(RYDBERG_SUM)
        % �ǽ����Ե���ʹ��Phelps���޽����Ե���
        % ������
        
        %%%%%%% ����Phelps
        data_dir={[path_list.cross_sections '\e_H2_Phelps_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% ׼������ԭʼ����
        idx_use=[]; % ���ù������
        E_head=1e-3; % ����������Χ
        E_end=1e4; 
        
        %%%%%%% ���Դ���
        % ʹ��IST
        data_dir={[path_list.cross_sections '\e_H2_Phelps_20201214\e_H2_IST_20201214']};
        Xsec_IST = ImportLXCat;
        Xsec_IST=Xsec_IST.init({'e-H2-IST'}, data_dir, 0);
        Xsec.ela=Xsec_IST.ela;
        Xsec.info.ela=Xsec_IST.info.ela;
        Xsec.eff{1}{1}=[0 0];
        
        %%%%%%% ��������
        idx_use=[1:4 6:10,12,14];
        Xsec=Xsec.choose_exc(1, idx_use);
        
        %%%%%%% ��������
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %��info�������ֵ��Ϣ
        
    case '1994Ness'
        data_dir={[path_list.root '\examples\example_ExB_drift']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        Xsec.tot{1,1}=Xsec.ela{1,1}{1,1}(:,2);
        Xsec.energy=Xsec.ela{1,1}{1,1}(:,1);
    otherwise
        error('[ERROR] No such type of Xsec_e_H2')
end
% ���
disp(['[INFO] Get Xsec : ' Xsec.name{1} ' ok']);
if flag_display
    Xsec.output_info();
end
end

%% aid function
% ����
% function data_dir=get_data_dir(name_part, path_list)
% % û��Ҫר��дһ������ͬ�ļ�����׺���Ӻ���
% name_part=strrep(name_part,'-','_');
% data_dir={[path_list.cross_sections '\e_H2_Morgan_20201214']};
% end