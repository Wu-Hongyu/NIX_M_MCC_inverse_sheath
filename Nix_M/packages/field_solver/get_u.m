function [ u ] = get_u( u0,A,b,type_solver)
% 求解 泊松方程-稀疏矩阵方程
switch type_solver
    case 'direct inverse'
        u=full(A\b);%左除解poisson方程，然后去稀疏
    otherwise
        tolerance= 1e-8;%容差
        max_num_iteration=500;%最大迭代次数
        switch type_solver
            case 'CG'
                [u,flag,reltive_residual_error,num_iteration] =pcg(A,b,tolerance,max_num_iteration,[],[],u0);
            case 'ICCG'
%                 C = rcond(full(A))
                L = ichol(A);
                [u,flag,reltive_residual_error,num_iteration] =pcg(A,b,tolerance,max_num_iteration,L,L',u0);
            case 'BiCGstab'
                % MATLAB官方AMG预处理器仅有并行版。FileExchange有一个MG代码
                % https://ww2.mathworks.cn/help/parallel-computing/solve-differential-equation-using-multigrid-preconditioner-on-distributed-discretization.html
                [u,flag,reltive_residual_error,num_iteration] =bicgstab(A,b,tolerance,max_num_iteration,[],[],u0);
            otherwise
                error('Not Done');
                % 并行求解器
                % https://ww2.mathworks.cn/help/parallel-computing/Use-Distributed-Arrays-to-Solve-Systems-of-Linear-Equations-with-Iterative-Methods.html
        end
        if (flag)
            error([type_solver '未收敛，flag=' num2str(flag)])
        else
            fprintf('%s %d次迭代后残差%d < 容差%d\n',type_solver,num_iteration,reltive_residual_error,tolerance)
        end
end
end