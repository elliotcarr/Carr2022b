function [mu,alpha] = weibull_model(d,Lbnd,Rbnd)
% Computes the Weibull model parameters mu and alpha [Sections 5.1-5.3 in manuscript]

a0 = Lbnd{1};
b0 = Lbnd{2};
L0 = Lbnd{3};
a1 = Rbnd{1};
b1 = Rbnd{2};
L1 = Rbnd{3};

% Pade coefficients
p = [0.45810,0.15757,1.49126,0.13963,-1.31348,3.28085];

if L0 > 0 % annular/spherical-shell system
    
    syms s r u w;
    
    I1 = integrals(1,d,L0,L1);
    I2 = integrals(2,d,L0,L1);
    I3 = integrals(3,d,L0,L1);
    I4 = integrals(4,d,L0,L1);
    eta = a0*(a1*I1+b1*L1^(1-d)) + a1*b0*L0^(1-d);
    beta(1) = ((a1*I1+b1*L1^(1-d))*(a0*L0^2-2*b0*L0) + b0*L0^(1-d)*(a1*L1^2+2*b1*L1)) / eta;
    beta(2) = (a0*a1*(L1^2-L0^2) + 2*(a0*b1*L1+a1*b0*L0)) / eta;
    gamma1s(1) = -a0*(L0^4/(4*(d+2)) - beta(1)*L0^2/(2*d)) + b0*(L0^3/(d+2) - beta(1)*L0/d);
    gamma1s(2) = -a1*(L1^4/(4*(d+2)) - beta(1)*L1^2/(2*d) - beta(2)*I3) - b1*(L1^3/(d+2)-beta(1)*L1/d-beta(2)*L1^(1-d)*I2);
    gamma1(1) = ((a1*I1+b1*L1^(1-d))*gamma1s(1) + b0*L0^(1-d)*gamma1s(2)) / eta;
    gamma1(2) = (a0*gamma1s(2)-a1*gamma1s(1)) / eta;    
    kappa = 2*(L1^d-L0^d)*((L1^(d+4)-L0^(d+4))/(4*(d+2)*(d+4)) - beta(1)*(L1^(d+2)-L0^(d+2))/(2*d*(d+2)) + ...
        gamma1(1)*(L1^d-L0^d)/d - beta(2)*I4 + gamma1(2)*I2) / ...
        ((beta(1)*(L1^d-L0^d)/d - (L1^(d+2)-L0^(d+2))/(d+2) + beta(2)*I2)^2); 
    alpha = (p(5)*kappa-p(2) - sqrt((p(5)*kappa-p(2))^2-4*(p(3)-p(6)*kappa)*(p(1)-p(4)*kappa))) / ...
        (2*(p(3)-p(6)*kappa));
    mu = alpha*(beta(1)*(d+2)*(L1^d-L0^d) + beta(2)*d*(d+2)*I2 - d*(L1^(d+2)-L0^(d+2))) / ...
        (2*d*(d+2)*(L1^d-L0^d)*gamma(1/alpha));
    
else % circular/spherical system
    
    sigma1 = b1/a1; L = L1;

    if sigma1 > 0 % semi-absorbing boundary
        
        kappa = (d+2)*(2*L^4 + sigma1*(d+4)*(2*L^3 + sigma1*(d+2)*L^2))/ ...
            ((d+4)*(L^2 + sigma1*(d+2)*L)^2);
        alpha = (p(5)*kappa-p(2) - sqrt((p(5)*kappa-p(2))^2-4*(p(3)-p(6)*kappa)*(p(1)-p(4)*kappa))) / ...
            (2*(p(3)-p(6)*kappa));
        mu = alpha*(L^2 + sigma1*(d+2)*L)/(d*(d+2)*gamma(1/alpha));
        
    else % absorbing boundary
        
        if d == 1
            alpha = 0.84883;
        elseif d == 2
            alpha = 0.78258;
        elseif d == 3
            alpha = 0.74510;
        end
        
        mu = alpha*L^2/(d*(d+2)*gamma(1/alpha)); 
        
    end

end