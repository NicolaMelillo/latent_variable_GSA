function Y = my_reg(X, rho, flg)

if rho ~=0
    lambda1 = sqrt(abs(rho));
    s_eps = sqrt(1-lambda1^2);
    eta = X(:,5);
else
    lambda1 = 0;
    %lambda2 = 0;
    eta = 0;
    s_eps = 1;
end

X1 = +lambda1*eta + s_eps*X(:, 1);
X2 = X(:, 2);
X3 = X(:, 3);

if flg == 1
    Y = X1 + X2 + X2.*X3;
elseif flg == 2
    Y = X1 + X2 + X1.*X3;
elseif flg == 3
    lambda2 = sign(rho) * sqrt(abs(rho));
    X4 = +lambda2*eta + s_eps*X(:, 4);
    Y = X1 + X2 + X3 + X4;
end
    
end

