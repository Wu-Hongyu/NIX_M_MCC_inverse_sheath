function flag=isequal_with_error(new_value,true_value,type,err)
% 判断两个数是否在给定容差内相等
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
            error('无法与0比较相对误差')
        else
            if abs((new_value-true_value)/true_value)<err
                flag=true;
            else
                flag=false;
            end
        end
end
end