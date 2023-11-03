function [Ps,t] = stochastic_model(d,P,P0,P1,delta,tau,L0,L1,Ns,Np,T)
% Implements the stochastic model [Section 2 in manuscript]

Nsteps = round(T/tau);

Ps = zeros(Ns,Nsteps+1);
Ps(:,1) = 1;
for k = 1:Ns % loop over stochastic simulations
    
    % initial uniform density
    r = nthroot(rand(Np,1)*(L1^d - L0^d) + L0^d,d);
    if d == 1
        xi = r;
    elseif d == 2
        theta = 2*pi*rand(Np,1);
        xi = [r.*cos(theta), r.*sin(theta)];
    elseif d == 3
        theta = 2*pi*rand(Np,1);
        phi = acos(2*rand(Np,1)-1);
        xi = [r.*cos(theta).*sin(phi), r.*sin(theta).*sin(phi), r.*cos(phi)];
    end
    
    x = xi;
    released = zeros(Np,1);

    for j = 1:Np % loop over particles
        
        for i = 1:Nsteps
        
            rn = rand;
            if rn <= P % move
                if d == 1
                    dx = delta*sign(rand-0.5);
                elseif d == 2
                    theta = 2*pi*rand;
                    dx = delta*[cos(theta),sin(theta)];
                elseif d == 3
                    theta = 2*pi*rand;
                    phi = acos(2*rand-1);
                    dx =  delta*[cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];
                end
                xt = x(j,:) + dx;
                r = norm(xt,2);
            else
                xt = x(j,:);
                r = norm(xt,2);
            end
            
            % inner boundary
            if r <= L0
                rn = rand;
                if rn <= P0 % absorb
                    x(j,:) = xt;
                    released(j) = i;
                    break;
                else
                    continue;
                end
            end
            
            % outer boundary
            if r >= L1
                rn = rand;
                if rn <= P1 % absorb
                    x(j,:) = xt;
                    released(j) = i;
                    break;
                else
                    continue;
                end
            end
            
            x(j,:) = xt;
            
        end
        
    end
    
    % calculate proportion of particles remaining in system
    for i = 1:Nsteps
        Nrp = sum(released<=i);
        Ps(k,i+1) = 1-Nrp/Np;
    end
    
    t = tau*[0:Nsteps];
    fprintf('%g %%\n',round(100*(k/Ns)))

end
