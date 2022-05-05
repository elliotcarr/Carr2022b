function I = simpsons_rule(h,f)
% Simpson's rule approximation to integral of f(x) from x = 0 to x = h*(length(f)-1)
n = length(f);
if mod(n,2) == 0
    error('n must be odd');
end
I = h/3*(f(1) + 4*sum(f(2:2:n-1)) + 2*sum(f(3:2:n-1)) + f(n));
