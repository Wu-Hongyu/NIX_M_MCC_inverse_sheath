function flag=isequal_with_error(new_value,true_value,type,err)
% �ж��������Ƿ��ڸ����ݲ������
switch type
    case 'AbsTol'
        if abs(new_value-true_value)<err
            flag=true;
        else
            flag=false;
        end
    case 'RelTol'
        if 0==true_value
            flag=false;
            error('�޷���0�Ƚ�������')
        else
            if abs((new_value-true_value)/true_value)<err
                flag=true;
            else
                flag=false;
            end
        end
end
end