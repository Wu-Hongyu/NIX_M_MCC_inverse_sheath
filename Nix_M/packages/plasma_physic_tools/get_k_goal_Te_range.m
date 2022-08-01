function [Te_array,k_array]=get_k_goal_Te_range(xsec)
% ���Ŀ��Te��Χ�ڵ�k_array
    Te_min=0.1;
    Te_max=100;
    num_k=1E2;
    Te_array=logspace(log10(Te_min),log10(Te_max),num_k);
    k_array=zeros(1,num_k);
    for j=1:num_k
        Te_j=Te_array(j);
        k_array(j)=get_k_e_heavy(xsec, Te_j);
    end
end

%PengChen: test ok 2020/12/24 17:13:44 
% ��LZS��e-H exc 1to2���ʶԱȣ�����һ��
% ��һ�¿�������Ϊ�Ե�����ֵ����Ĵ���ͬ