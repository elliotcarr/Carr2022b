function I = integrals(n,d,L0,L1)

if n == 1
    if d == 1
        I = L1 - L0;
    elseif d == 2
        I = log(L1/L0);
    elseif d == 3
        I = (L1-L0)/(L0*L1);
    end
elseif n == 2
    if d == 1
        I = (L1-L0)^2/2;
    elseif d == 2
        I = (2*L1^2*log(L1/L0)-(L1^2-L0^2))/4;
    elseif d == 3
        I = (L1-L0)^2*(L0+2*L1)/(6*L0);
    end
elseif n == 3
    if d == 1
        I = (L1-L0)^3/6;
    elseif d == 2
        I = ((L0^2+L1^2)*log(L1/L0)-(L1^2-L0^2))/4;
    elseif d == 3
        I = (L1-L0)^3/(6*L0*L1);
    end
elseif n == 4
    if d == 1
        I = (L1-L0)^4/24;
    elseif d == 2
        I = (4*L1^2*(L1^2+2*L0^2)*log(L1/L0)+L0^4-5*L1^4+4*L0^2*L1^2)/64;
    elseif d == 3
        I = (L1-L0)^4*(L0+4*L1)/(120*L0);
    end
end