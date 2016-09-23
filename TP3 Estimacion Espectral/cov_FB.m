function R = cov_FB (y, L)

    N = length(y);
    R_hat = zeros(L,L);  %inicializacion
    
    for i = L : N
        R_hat = R_hat + y(i : -1 : i-L+1) * (y(i : -1 : i-L+1)');
    end
    
    R_hat = R_hat / N ;
    
    J = eye(L);
    J = rot90(J);
    
    R = (R_hat + J * R_hat.' * J) / 2;
    
end