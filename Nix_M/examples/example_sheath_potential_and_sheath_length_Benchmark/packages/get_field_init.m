function [ u, A, b_extra ] = get_field_init( simulation )
n=simulation.num_grid_point;
% �����λ�߽���������װϵ������
u=zeros(n,1);%��λ��ֵ
leading_diagonal=-2*ones(n,1);
lower_diagonal=ones(n,1);
upper_diagonal=ones(n,1);
% spdiags��װn �� n ��ϡ��ϵ������
% λ�����Խ����·��ĶԽ������ȴ��еĶ�����ȡԪ�ء�
% λ�����Խ����Ϸ��ĶԽ������ȴ��еĵײ���ȡԪ�ء�
b_extra=zeros(2,1);%���Ҷ������Ķ������
switch simulation.field_boundaries_type
    case 0 %00����ĩ��һ��߽�
        u(1)=simulation.field_boundaries(1);
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        u(end)=simulation.field_boundaries(2);
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    case 1 %01���׵�һ�࣬ĩ�ڶ���
        u(1)=simulation.field_boundaries(1);
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        lower_diagonal(end-1)=2;
        b_extra(end)=-2*simulation.field_boundaries(2)*simulation.dx;
    case 2 %10���׵ڶ��࣬ĩ��һ��
        upper_diagonal(2)=2;
        b_extra(1)=2*simulation.field_boundaries(1)*simulation.dx;
        
        u(end)=simulation.field_boundaries(2);
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    case 3 %11�����ڱ߽�����
        disp('ʹ�õ糡���ڱ߽��������߽��λȡΪ��ο���λ������simulation.field_boundaries')
        u(1)=0;
        leading_diagonal(1)=1;
        upper_diagonal(2)=0;
        b_extra(1)=u(1);
        
        u(end)=0;
        leading_diagonal(end)=1;
        lower_diagonal(end-1)=0;
        b_extra(end)=u(end);
    otherwise
        error('Not Done');
end
A=spdiags([lower_diagonal,leading_diagonal,upper_diagonal],[-1,0,1],n,n);
% ��С�Ҷ�����������
constants=get_constants();% ȫ�ֳ��� �ṹ��
coeff_temp=-constants.eps0/(simulation.dx*simulation.dx);
A=coeff_temp*A;
b_extra=coeff_temp*b_extra;
end
