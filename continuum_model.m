function [t,c] = continuum_model(D,Nr,d,Lbnd,Rbnd,T,Nt)
% Solves the continuum model [Equations (2)-(5) in manuscript] using finite volume method

a0 = Lbnd{1};
b0 = Lbnd{2};
L0 = Lbnd{3};
a1 = Rbnd{1};
b1 = Rbnd{2};
L1 = Rbnd{3};

dt = T/Nt; % time step size
t = linspace(0,T,Nt+1)'; % discrete times
h = (L1-L0)/(Nr-1); % node spacing
r = linspace(L0,L1,Nr); % nodes
rw(1) = r(1); rw(2:Nr) = (r(1:Nr-1)+r(2:Nr))/2; % west boundaries 
re(1:Nr-1) = rw(2:Nr); re(Nr) = r(Nr); % east boundaries

% Finite volume space discretisation
A = zeros(Nr,Nr); b = zeros(Nr,1);
if b0 ~= 0 % reflecting or semi-absorbing boundary
    A(1,1) = -(D*re(1)^(d-1)/h + D*rw(1)^(d-1)*a0/b0);
    A(1,2) = D*re(1)^(d-1)/h;
    if L0 > 0 || d == 1
        A(1,1) = A(1,1)/(h/2*r(1)^(d-1));
        A(1,2) = A(1,2)/(h/2*r(1)^(d-1));
        b(1) = b(1)/(h/2*r(1)^(d-1));
    end
end
for i = 2:Nr-1 % interior nodes
    A(i,i-1) = D*rw(i)^(d-1)/h;
    A(i,i) = -(D*re(i)^(d-1) + D*rw(i)^(d-1))/h;
    A(i,i+1) = D*re(i)^(d-1)/h;
    A(i,i-1) = A(i,i-1)/(h*r(i)^(d-1));
    A(i,i) = A(i,i)/(h*r(i)^(d-1));
    A(i,i+1) = A(i,i+1)/(h*r(i)^(d-1));
end
if b1 ~= 0 % reflecting or semi-absorbing boundary
    A(Nr,Nr-1) = D*rw(Nr)^(d-1)/h;
    A(Nr,Nr) = -(D*rw(Nr)^(d-1)/h + D*re(Nr)^(d-1)*a1/b1);
    A(Nr,Nr-1) = A(Nr,Nr-1)/(h/2*r(Nr)^(d-1));
    A(Nr,Nr) = A(Nr,Nr)/(h/2*r(Nr)^(d-1));
    b(Nr) = b(Nr)/(h/2*r(Nr)^(d-1));
end

% Crank Nicolson time descretisation
theta = 1/2;
if b0 == 0 % absorbing boundary
    At(1,:) = [1,zeros(1,Nr-1)];
else
    At(1,:) = -A(1,:);
    if L0 > 0 || d == 1
        At(1,:) = [1,zeros(1,Nr-1)]-theta*dt*A(1,:);
    end
end
for i = 2:Nr-1 % interior nodes
    At(i,i-1) = -theta*dt*A(i,i-1);
    At(i,i) = 1-theta*dt*A(i,i);
    At(i,i+1) = -theta*dt*A(i,i+1);
end
if b1 == 0 % absorbing boundary
    At(Nr,:) = [zeros(1,Nr-1),1];
else
    At(Nr,:) = [zeros(1,Nr-1),1]-theta*dt*A(Nr,:);
end
At = sparse(At);

c = zeros(Nr,Nt+1); 
c(:,1) = ones(Nr,1); % initial uniform concentration

% Time stepping
for n = 1:Nt
    bt = c(:,n) + (1-theta)*dt*A*c(:,n) + dt*b;
    if b0 == 0 % absorbing boundary
        bt(1) = 0;
    else
        if (L0 == 0 && d ~= 1)
            bt(1) = b(1);
        end
    end
    if b1 == 0
        bt(Nr) = 0;
    end
    c(:,n+1) = At\bt;
end
