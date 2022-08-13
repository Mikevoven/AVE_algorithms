function [x_k,i,err] = MBB(A,B,b)
%MBB 方法用于求解广义绝对值方程：Ax+B|x|=b
%   在运行该文件之前，建议先运行 MBBCheck 检查是否满足 MBB 方法的使用条件

%   Version 1.0
%   Copyright 2022 S. Yang

tol = 1e-6;
n = 500;

x_k_1 = zeros(length(b),1);
inv_alpha_0 = 1;
beta = 0.78;
sigma = 0.001;
m = 1;

z_0 = x_k_1 - beta^m*inv_alpha_0*(A*x_k_1 + B*abs(x_k_1) - b);

g_z0 = A*z_0 + B*abs(z_0) - b;
zeta_0 = (g_z0' *(x_k_1 - z_0))/(g_z0' * g_z0);
x_k = x_k_1 - zeta_0*g_z0;

for i = 1:n
    err = norm(A*x_k + B*abs(x_k) - b);
    if err > tol
        g_k_1 = A*x_k_1 + B*abs(x_k_1) - b;
        g_k = A*x_k + B*abs(x_k) - b;
        inv_alpha_k = (x_k - x_k_1)'*(x_k - x_k_1)/((x_k - x_k_1)'*(g_k - g_k_1));
        
        % update m
        m = 0;
        try_point = x_k - beta^m*inv_alpha_k*g_k;
        g_try_point = A*try_point + B*abs(try_point) -b;
        while g_try_point'*g_k < sigma*beta^m*(g_k'*g_k)*inv_alpha_k
            m = m+1;
            try_point = x_k - beta^m*inv_alpha_k*g_k;
            g_try_point = A*try_point + B*abs(try_point) -b;
        end
        
        % update x_{k+1}
        z_k = x_k - beta^m*inv_alpha_k*g_k;
        g_zk = A*z_k + B*abs(z_k) - b;
        zeta_1 =  (g_zk' *(x_k - z_k))/(g_zk' * g_zk);

        [x_k_1, x_k] = deal(x_k, x_k - zeta_1*g_zk);
    else
        break
    end
end
end

