function [x,S,err,iter,flag] = SignAccord(A, B, b)
%Finds a solution to Ax + B|x| = b or states
%   singularity of [A − |B|, A + |B|].
x = [];
S = [];
flag = 'singular';
if det(A) == 0
    S = A;
    return
end

z = sign(A\b);
if det(A+B*diag(z)) == 0
    S = A+B*diag(z);
    return
end

n = length(b);
p = sparse(n,1);
x = (A+B*diag(z))\b;
C = -(A+B*diag(z))\B;
iter = 0;
while any(x.*z < 0)
    iter = iter+1;
    k = find(x.*z < 0, 1);
    E_k = sparse(k,k,1,n,n);
    if 1+2*z(k)*C(k,k) <= 0
        S = A + B*(diag(z)+1/(C(k,k))*E_k);
        x = [];
        return
    end
    p = p+1;
    if log(p)/log(2) > n-k
        x = [];
        return
    end
    z(k) = -z(k);
    alpha = 2*z(k)/(1 - 2*z(k)*C(k,k));
    x = x + alpha*x(k)*C(:,k);
    C = C+alpha*C(:,k)*C(k,:);
end
flag = 'solution';
err = norm(A*x+B*abs(x)-b);