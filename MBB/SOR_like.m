function [x0,i,err] = SOR_like(A,B,b,w)
tol = 10^(-6);               % 误差界
N = 500;                    % 最大迭代次数

x0 = zeros(length(b),1);
y0 = B*abs(x0);
for i = 1:N+1
    x0 = (1-w)*x0 + w*(A\(b-y0));
    y0 = (1-w)*y0 + w*B*abs(x0);
    err = norm(A*x0 + B*abs(x0) - b);
    if  err <= tol
        break
    elseif i == N
        disp('超出最大迭代次数')
    end
end
end
