function [Xsec] = get_Xsec_Ion(type, flag_display)
%GET_XSEC_HP_H2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
switch type
    case 'H+_H2'
        path_list=get_path_init('test000');
        recation_type={['Hp-H2' type]};  %�����cell
        data_dir={[path_list.cross_sections '\Hp_H2_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H+_H'
        path_list=get_path_init('test000');
        recation_type={['Hp-H' type]};  %�����cell
        data_dir={[path_list.cross_sections '\Hp_H_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H2+_H2'
        path_list=get_path_init('test000');
        recation_type={['H2p-H2' type]};  %�����cell
        data_dir={[path_list.cross_sections '\H2p_H2_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H2+_H'
        path_list=get_path_init('test000');
        recation_type={['H2p-H' type]};  %�����cell
        data_dir={[path_list.cross_sections '\H2p_H_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H3+_H2'
        path_list=get_path_init('test000');
        recation_type={['H3p-H2' type]};  %�����cell
        data_dir={[path_list.cross_sections '\H3p_H2_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H3+_H'
        path_list=get_path_init('test000');
        recation_type={['H3p-H' type]};  %�����cell
        data_dir={[path_list.cross_sections '\H3p_H_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H-_H2'
        path_list=get_path_init('test000');
        recation_type={['Hn-H2' type]};  %�����cell
        data_dir={[path_list.cross_sections '\Hn_H2_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
    case 'H-_H'
        path_list=get_path_init('test000');
        recation_type={['Hn-H' type]};  %�����cell
        data_dir={[path_list.cross_sections '\Hn_H_Xsec']};
        Xsec = ImportLXCat;
        Xsec=Xsec.init(recation_type, data_dir, flag_display);
        
        
end

