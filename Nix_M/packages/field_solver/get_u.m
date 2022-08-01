function [ u ] = get_u( u0,A,b,type_solver)
% ��� ���ɷ���-ϡ����󷽳�
switch type_solver
    case 'direct inverse'
        u=full(A\b);%�����poisson���̣�Ȼ��ȥϡ��
    otherwise
        tolerance= 1e-8;%�ݲ�
        max_num_iteration=500;%����������
        switch type_solver
            case 'CG'
                [u,flag,reltive_residual_error,num_iteration] =pcg(A,b,tolerance,max_num_iteration,[],[],u0);
            case 'ICCG'
%                 C = rcond(full(A))
                L = ichol(A);
                [u,flag,reltive_residual_error,num_iteration] =pcg(A,b,tolerance,max_num_iteration,L,L',u0);
            case 'BiCGstab'
                % MATLAB�ٷ�AMGԤ���������в��а档FileExchange��һ��MG����
                % https://ww2.mathworks.cn/help/parallel-computing/solve-differential-equation-using-multigrid-preconditioner-on-distributed-discretization.html
                [u,flag,reltive_residual_error,num_iteration] =bicgstab(A,b,tolerance,max_num_iteration,[],[],u0);
            otherwise
                error('Not Done');
                % ���������
                % https://ww2.mathworks.cn/help/parallel-computing/Use-Distributed-Arrays-to-Solve-Systems-of-Linear-Equations-with-Iterative-Methods.html
        end
        if (flag)
            error([type_solver 'δ������flag=' num2str(flag)])
        else
            fprintf('%s %d�ε�����в�%d < �ݲ�%d\n',type_solver,num_iteration,reltive_residual_error,tolerance)
        end
end
end