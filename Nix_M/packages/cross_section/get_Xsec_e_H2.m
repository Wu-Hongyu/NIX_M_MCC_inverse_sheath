function [ Xsec ] = get_Xsec_e_H2( type, flag_display)
% 处理ImportLXCat得到的e-H2原始数据，获取MCC所需截面
% 使用
% type：e_H2散射截面数据种类名，一般即数据库名
% flag_display：0则不画图不输出反应种类，1则显示
% 注意，本文件根据.\others\choose_Xsec_e_n_of_hydrogen_plasma
% \compare_Xsec_e_H2.m实现
% 物理性的比较在compare_Xsec_e_H2
% 功能性的test在test_get_Xsec

% 通过全局替换E_head和E_end值可以修改截面可用区间
path_list=get_path_init('test000');
recation_type={['e-H2-' type]};  %处理成cell
switch type       
    case '2020LZS' % 李增山
        % 弹性使用CCC
        % 忽略旋转激发
        % 振动激发使用2003Janev，选取v1,v2
        % 电子激发使用CCC，根据过程的N和截面选用13个过程，
        % LZS根据Fantz2006给定阈值
        % 电离使用2011Wuenderlich，选取非解离性电离。可扩展至解离性
        % 吸附使用2003Janev，选取v=0。可扩展至v=0~14，但负Eth被处理成E_head
        
        %%%%%%% 导入CCC原始数据
        data_dir={[path_list.cross_sections '\e_H2_CCC_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% 准备补充原始数据
        idx_use=[]; % 可用过程序号
        E_head=1e-3; % 可用能量范围
        E_end=1e4; 
        
        %%%%%%% 弹性处理
        for E_end_i=logspace(log10(400),log10(E_end),10)
            Xsec=Xsec.add_key_points2(E_head, E_end_i);
        end
        % CCC数据只到300eV，若只在E_end处补一个点，trapz时可能忽略掉了
        
        %%%%%%% 振动激发处理
        % 排在exc末尾，在选择时提前
        num0_exc=length(Xsec.exc{1});
        % 截面拟合表达式 2003Janev-ch 4.1.1-eq.74
        xsec_v=@(ksv,dE,E) 5.78*ksv/(dE^4)*((dE/E)^2)...
            *(1-dE/E)^(6.11/sqrt(dE))*1E-16*1E-4;
        % dE: ΔE，the energy difference of v = 0 and v'
        % ksv: factor
        % v=0->v=1
        Eth_v1=0.516; % 以ΔE为阈值
        xsec_v1=@(E) xsec_v(1,Eth_v1,E);
        x_v1=logspace(log10(Eth_v1),log10(E_end),100); 
        y_v1=arrayfun(xsec_v1,x_v1);
        % 更新Xsec
        k=num0_exc+1;
        Xsec.exc{1}{k}=[x_v1',y_v1'];
        Xsec.excThresh{1}{k}=Eth_v1;
        Xsec.info.exc{1}{k}='E+H2->E+H2(v=1),Excitation';
        idx_use=[idx_use k];
        % v=0->v=2
        Eth_v2=1.003; % 以ΔE为阈值
        xsec_v2=@(E) xsec_v(0.628,Eth_v2,E);
        x_v2=logspace(log10(Eth_v2),log10(E_end),100); 
        y_v2=arrayfun(xsec_v2,x_v2);
        % 更新Xsec
        k=num0_exc+2;
        Xsec.exc{1}{k}=[x_v2',y_v2'];
        Xsec.excThresh{1}{k}=Eth_v2;
        Xsec.info.exc{1}{k}='E+H2->E+H2(v=2),Excitation';
        idx_use=[idx_use k];
        
        %%%%%%% 电子激发处理
        % CCC exc 1为GTSC，总截面
        % 第一类ele：LZS根据截面与k大小选取，为CCC部分ele
        % LZS以2006Fantzb电子态能量作为阈值
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
        % 分解
        % 电子撞击导致的氢分子分解反应，一般经由激发，
        % 即先发生激发反应，然后再激发态氢分子分解
        % 由于MCC不模拟H2/H，因此可以不考虑第二步
        % 此处补充考虑b3电子激发，1987Janev的阈值
        Xsec.excThresh{1}{20}=8.5; %b3
        idx_use=[idx_use 20];

        % 第二类ele：PengChen根据k*ΔE选取，即CCC全部ele
%         % 根据2003Janev给定阈值，其中h3和J1使用2006Fantzb电子态能量
%         p_array=[2:5 7:8 10:20];
%         Eth_array=[11.72, 15.47, 14.85, 12.754, 13.29, 11.72, 15.555, 14.996, 13.6, 13, 13.13, 14.816, 14.98, 13.98, 14.824, 14.07, 7.93];
%         for j=1:length(p_array)
%             Xsec.excThresh{1}{p_array(j)}=Eth_array(j);
%         end
%         idx_use=[idx_use p_array];

        % 选择反应
        Xsec=Xsec.choose_exc(1, idx_use);
        
        for E_end_i=logspace(log10(300),log10(E_end),10)
            Xsec=Xsec.add_key_points_exc(E_head, E_end_i);
        end
        
        %%%%%%% 电离处理
        % 据LZS，参考2011Wuenderlich
        % 应该是2011Wuenderlich理论计算截面后，LZS对Emin到1e3eV做拟合
        % 参数
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
        Emin_ion=[16.01 18.47 41.98]; % 对应于非0截面 的Emin？
        Emax_ion=1e3;
        xsec_ion=@(E,i) sum( coeff_ion(:,i).*(log(E)).^(10:-1:0)' );
        % 电离-非解离性
        Eth_ion1=15.4; % LZS文档上的阈值
        xsec_ion1=@(E) xsec_ion(E,1);
        x_ion1=logspace(log10(Emin_ion(1)),log10(Emax_ion),100); 
        y_ion1=arrayfun(xsec_ion1,x_ion1);
        % 更新Xsec
        k=1;
        Xsec.ion{1}{k}=[x_ion1',y_ion1'];
        Xsec.ionThresh{1}{k}=Eth_ion1;
        Xsec.info.ion{1}{k}='E+H2->2E+H2+,Ionization';
        % 电离-解离性 g 截面小于1e22，且阈值较高，忽略
        Eth_ion2=18.15; % LZS文档上的阈值
        xsec_ion2=@(E) xsec_ion(E,2);
        x_ion2=logspace(log10(Emin_ion(2)),log10(Emax_ion),100); 
        y_ion2=arrayfun(xsec_ion2,x_ion2);
        % 更新Xsec
        k=k+1;
        Xsec.ion{1}{k}=[x_ion2',y_ion2'];
        Xsec.ionThresh{1}{k}=Eth_ion2;
        Xsec.info.ion{1}{k}='E+H2->2e+H2+(2Σg+)->2E+H+H+,Ionization';
        % 电离-解离性 u
        Eth_ion3=30.6; % LZS文档上的阈值
        xsec_ion3=@(E) xsec_ion(E,3);
        x_ion3=logspace(log10(Emin_ion(3)),log10(Emax_ion),100); 
        y_ion3=arrayfun(xsec_ion3,x_ion3);
        % 更新Xsec
        k=k+1;
        Xsec.ion{1}{k}=[x_ion3',y_ion3'];
        Xsec.ionThresh{1}{k}=Eth_ion3;
        Xsec.info.ion{1}{k}='E+H2->2e+H2+(2Σu+)->2E+H+H+,Ionization';
        % 通过与其他数据源比较，test ok
        for E_end_i=logspace(log10(Emax_ion),log10(E_end),10)
            Xsec=Xsec.add_key_points_ion(E_head, E_end_i);
        end
        
        %%%%%%% 吸附处理
        % 系数 2003Janev-Tab 27
        Eth_att=[3.72, 3.21, 2.72, 2.26, 1.83, 1.43, 1.36, 0.713, 0.397, 0.113, -0.139, -0.354, -0.529, -0.659, -0.736];
        sigma_max_att=[0.0000322, 0.000518, 0.00416, 0.022, 0.122, 0.453, 1.51, 4.48, 10.1, 13.9, 11.8, 8.87, 7.11, 5, 3.35];
        E0_att=0.45;
        % 截面拟合表达式 2003Janev-ch 4.4.1-eq.124
        xsec_att=@(E,v) sigma_max_att(v+1)*1E-20*exp(-(E-abs(Eth_att(v+1)))/E0_att);
%         for v=0:14
        for v=0:0 %暂时只考虑v=0振动态
            xsec_att_v=@(E) xsec_att(E,v);
            Eth_att_v=Eth_att(v+1);
            if Eth_att_v>0
                x_att_v=logspace(log10(Eth_att_v+1e-3),log10(E_end),100);
            else
                Eth_att_v=E_head;
                x_att_v=logspace(log10(E_head),log10(E_end),100);
            end
            y_att_v=arrayfun(xsec_att_v,x_att_v);
            % 更新Xsec
            k=v+1;
            Xsec.att{1}{k}=[x_att_v',y_att_v'];
            Xsec.attThresh{1}{k}=Eth_att_v;
            Xsec.info.att{1}{k}=['E+H2(v=' num2str(v) ')->H+H-,Attachment'];
        end
        % test结果见Nix_M截面与文献对比.pptx
        
        %%%%%%% 其他处理
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %在info中添加阈值信息
    case 'Phelps-m'
        % Phelps各过程的Ref细节见pengchen笔记
        % 弹性使用IST（Phelps的ela由eff换算，有bug）
        % 旋转激发使用Phelps，有j0-2，j1-3
        % 振动激发使用Phelps，选用v1,2,3
        % 电子激发使用Phelps，选用b3,B1,c3,a3,C1,D3,N≥4(RYDBERG_SUM)
        % 非解离性电离使用Phelps，无解离性电离
        % 无吸附
        
        %%%%%%% 导入Phelps
        data_dir={[path_list.cross_sections '\e_H2_Phelps_20201214']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        %%% 准备补充原始数据
        idx_use=[]; % 可用过程序号
        E_head=1e-3; % 可用能量范围
        E_end=1e4; 
        
        %%%%%%% 弹性处理
        % 使用IST
        data_dir={[path_list.cross_sections '\e_H2_Phelps_20201214\e_H2_IST_20201214']};
        Xsec_IST = ImportLXCat;
        Xsec_IST=Xsec_IST.init({'e-H2-IST'}, data_dir, 0);
        Xsec.ela=Xsec_IST.ela;
        Xsec.info.ela=Xsec_IST.info.ela;
        Xsec.eff{1}{1}=[0 0];
        
        %%%%%%% 激发处理
        idx_use=[1:4 6:10,12,14];
        Xsec=Xsec.choose_exc(1, idx_use);
        
        %%%%%%% 其他处理
        Xsec=Xsec.init_part2(E_head, E_end);
        Xsec=Xsec.add_threshold_to_info(); %在info中添加阈值信息
        
    case '1994Ness'
        data_dir={[path_list.root '\examples\example_ExB_drift']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init_part1(recation_type, data_dir, flag_display);
        Xsec.tot{1,1}=Xsec.ela{1,1}{1,1}(:,2);
        Xsec.energy=Xsec.ela{1,1}{1,1}(:,1);
    otherwise
        error('[ERROR] No such type of Xsec_e_H2')
end
% 输出
disp(['[INFO] Get Xsec : ' Xsec.name{1} ' ok']);
if flag_display
    Xsec.output_info();
end
end

%% aid function
% 搁置
% function data_dir=get_data_dir(name_part, path_list)
% % 没必要专门写一个允许不同文件名后缀的子函数
% name_part=strrep(name_part,'-','_');
% data_dir={[path_list.cross_sections '\e_H2_Morgan_20201214']};
% end