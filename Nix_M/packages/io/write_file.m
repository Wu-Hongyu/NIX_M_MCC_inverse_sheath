function write_file( filename, type, output )
% д�ļ�
% type: �ο� fopen-������� permission - �ļ���������
switch type
    case 'char2wt' % д�ַ����鵽�ı��ļ�
        out_file=fopen(filename,'wt');
        if out_file<0
            error(['[ERROR] can not write ' filename])
        end
        fprintf(out_file, '%s', output);
        % ע�⣬output�е�\r\n���޷�ʶ���
        fclose(out_file);
    otherwise
        error('[ERROR] not done')
end
end