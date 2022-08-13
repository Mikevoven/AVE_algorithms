function [min_iter,opt_w] = SOR_like_w(A,B,b)
tol = 10^(-6);               % 误差界
N = 500;                    % 最大迭代次数

cnt = 0;
Iters = zeros(1,19);

for w = 0.1:.1:1.9
    cnt = cnt + 1;
    x0 = zeros(length(b),1);
    y0 = B*abs(x0);
    for i = 1:N
        x0 = (1-w)*x0 + w*(A\(b-y0));
        y0 = (1-w)*y0 + w*B*abs(x0);
        err = norm(A*x0 + B*abs(x0) - b);
        if  err <= tol
            break
        end
    end
    Iters(cnt) = i; 
end
    [min_iter,opt_w] = min(Iters);
    opt_w = opt_w/10;
    plot(Iters);
end
