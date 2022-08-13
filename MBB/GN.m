function [x_k,i,err] = GN(A,B,b)
tol = 10^(-6);               % 误差界
N = 500;                    % 最大迭代次数

x_k = zeros(length(b),1);
for i = 1:N+1
    x_k = (A + B*diag(sign(x_k)))\b;
    err = norm(A*x_k + B*abs(x_k) - b);
    if  err <= tol
        break
    elseif i == N
        disp('超出最大迭代次数')
    end
end
end

