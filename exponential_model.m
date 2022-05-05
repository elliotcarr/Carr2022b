function lambda = exponential_model(d,Lbnd,Rbnd)
% Computes the exponential model parameter lambda [Sections 4.1-4.3 in manuscript]

a0 = Lbnd{1};
b0 = Lbnd{2};
L0 = Lbnd{3};
a1 = Rbnd{1};
b1 = Rbnd{2};
L1 = Rbnd{3};

if L0 > 0 % annular/spherical-shell system
    
    syms s r; 
    I1 = integrals(1,d,L0,L1);
    I2 = integrals(2,d,L0,L1);
    eta = a0*(a1*I1+b1*L1^(1-d)) + a1*b0*L0^(1-d);
    beta(1) = ((a1*I1+b1*L1^(1-d))*(a0*L0^2-2*b0*L0) + b0*L0^(1-d)*(a1*L1^2+2*b1*L1)) / eta;
    beta(2) = (a0*a1*(L1^2-L0^2) + 2*(a0*b1*L1+a1*b0*L0)) / eta;
    lambda = (beta(1)*(d+2)*(L1^d-L0^d) + beta(2)*d*(d+2)*I2 - d*(L1^(d+2)-L0^(d+2))) / ...
        (2*d*(d+2)*(L1^d-L0^d));
    
else % circular/spherical system
    
    sigma1 = b1; L = L1;
    if sigma1 > 0 % semi-aborbing boundary
        lambda = (L^2 + sigma1*(d+2)*L)/(d*(d+2));
    else % absorbing boundary
        lambda = L^2/(d*(d+2));
    end

end