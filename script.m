close all
clear all
clc

%%
addpath('/Users/carre/Dropbox/Documents/Research/Code/Other packages/export_fig-master')
path_name = '/Users/carre/Dropbox/Documents/Research/Projects/Spatial averaged model/Paper PA/Figures/';

% Select test case
% Case = 'A'; % disc/sphere with absorbing boundary
% Case = 'B'; % disc/sphere with semi-absorbing boundary
% Case = 'C'; % annulus/spherical shell with reflecting inner and absorbing outer boundary
% Case = 'D'; % annulus/spherical shell with reflecting inner and semi-absorbing outer boundary
% Case = 'E'; % annulus/spherical shell with absorbing inner and outer boundary
Case = 'F'; % annulus/spherical shell with semi-absorbing inner and outer boundary

% save_figs = false;
save_figs = true;

stochastic = true;
% stochastic = false;

% Parameters
Nr = 501; % number of nodes
Nt = 1e4; % number of time steps
Ns1 = 100; Ns2 = 100; % number of stochastic random walk simulations
Np1 = 50; Np2 = 500; % number of particles
k = 2; % end time T satisfies Pw(T) = 10^(-k)

% Plot formatting
font_size = 22;
line_width = 4;
line_width_axis = 1.0;
colors = [0, 47, 108; 200, 16, 46; 255, 199, 44]/255;
colors(4,:) = (1-0.2*(1-colors(1,:)));
colors(5,:) = (1-0.4*(1-colors(1,:)));
font_interpreter = 'LaTeX';

erre = zeros(1,3); errw = zeros(1,3);
for d = 2:3 % dimensions
    
    P = 1; % probability of moving
    delta = 1; % step distance
    tau = 1; % step duration
    
    if isequal(Case,'A')
        P0 = 0; % probability of absorbing at inner boundary
        P1 = 1; % probability of absorbing at outer boundary
        L0 = 0; L1 = 100; % inner and outer boundary locations
    elseif isequal(Case,'B')
        P0 = 0; % probability of absorbing at inner boundary
        P1 = 0.2; % probability of absorbing at outer boundary
        L0 = 0; L1 = 100; % inner and outer boundary locations
    elseif isequal(Case,'C')
        P0 = 0; % probability of absorbing at inner boundary
        P1 = 1; % probability of absorbing at outer boundary
        L0 = 50; L1 = 100; % inner and outer boundary locations
    elseif isequal(Case,'D')
        P0 = 0; % probability of absorbing at inner boundary
        P1 = 0.2; % probability of absorbing at outer boundary
        L0 = 50; L1 = 100; % inner and outer boundary locations
    elseif isequal(Case,'E')
        P0 = 1; % probability of absorbing at inner boundary
        P1 = 1; % probability of absorbing at outer boundary
        L0 = 50; L1 = 100; % inner and outer boundary locations
    elseif isequal(Case,'F')
        P0 = 0.5; % probability of absorbing at inner boundary
        P1 = 0.2; % probability of absorbing at outer boundary
        L0 = 50; L1 = 100; % inner and outer boundary locations
    end
    
    % Continuum model parameters
    D = P*delta^2/(2*d*tau);
    if P0 == 1 % absorbing
        a0 = 1; b0 = 0;
    elseif P0 == 0 % reflecting
        a0 = 0; b0 = 1;
    else % semi-absorbing
        a0 = 1; b0 = delta/P0;
    end
    if P1 == 1 % absorbing
        a1 = 1; b1 = 0;
    elseif P1 == 0 % reflecting
        a1 = 0; b1 = 1;
    else % semi-absorbing
        a1 = 1; b1 = delta/P1;
    end
    Lbnd = {a0,b0,L0};
    Rbnd = {a1,b1,L1};
    
    %% Exponential model parameters
    lambda = exponential_model(d,Lbnd,Rbnd);
    
    %% Weibull model parameters
    [mu,alpha] = weibull_model(d,Lbnd,Rbnd);
    
    %% Continuum model
    %T = -lambda*log(tol)/D; % end time [T satisfies Pe(T) = 10^(-k)]
    T = mu*(k*log(10))^(1/alpha)/D; % end time [T satisfies Pw(T) = 10^(-k)]
    r = linspace(L0,L1,Nr)'; % node locations
    [t,c] = continuum_model(D,Nr,d,Lbnd,Rbnd,T,Nt);
    
    % Simpson's rule approximation to spatial average
    Pc = zeros(length(t),1);
    h = (L1-L0)/(Nr-1);
    for i = 1:length(t)
        Pc(i) = d*simpsons_rule(h,(r.^(d-1)).*c(:,i))/(L1^d-L0^d);
    end
    Pc(1) = 1;
    
    %% Exponential model
    Pe = exp(-D*t/lambda);
    
    %% Weibull model
    Pw = exp(-(D*t/mu).^alpha);
    
    %% Errors
    erre(d) = mean(abs(Pe-Pc)); % exponential
    errw(d) = mean(abs(Pw-Pc)); % weibull
    
    %% Stochastic simulations
    if stochastic
        [Ps1,ts1] = stochastic_model(d,P,P0,P1,delta,tau,L0,L1,Ns1,Np1,T);
        lower1 = quantile(Ps1,0.025); upper1 = quantile(Ps1,0.975);
        [Ps2,ts2] = stochastic_model(d,P,P0,P1,delta,tau,L0,L1,Ns2,Np2,T);
        lower2 = quantile(Ps2,0.025); upper2 = quantile(Ps2,0.975);
    end
    
    %% Plots
    figure;
    set(gcf,'Position',[434   239   560*0.8   520*0.7]);
    set(gcf,'Color','w')
    set(gcf,'Renderer','Painters');
    hold on
    if stochastic
        patch([ts1,fliplr(ts1)]',[lower1,fliplr(upper1)]',1,'FaceColor',colors(4,:),'EdgeColor',colors(4,:))
        patch([ts2,fliplr(ts2)]',[lower2,fliplr(upper2)]',1,'FaceColor',colors(5,:),'EdgeColor',colors(5,:))
    end
    plot(t,Pc,'-','LineWidth',line_width,'Color',colors(1,:));
    plot(t,Pe,'-','LineWidth',line_width,'Color',colors(2,:));
    plot(t,Pw,'--','LineWidth',line_width,'Color',colors(3,:));
    if isequal(font_interpreter,'LaTeX')
        xl = xlabel('$t$','Interpreter',font_interpreter);
    else
        xl = xlabel('t','Interpreter',font_interpreter);
    end
    set(xl,'Position',[T/2,-0.05])
    if isequal(font_interpreter,'LaTeX')
        yl = ylabel('$\mathcal{P}(t)$','Interpreter',font_interpreter);
    else
        yl = ylabel('P(t)','Interpreter',font_interpreter);
    end
    posyl = get(yl,'Position');
    set(yl,'Position',[posyl(1)+0.5*(T/50),0.5])
    xlim([0,T])
    ylim([0,1])
    leg.ItemTokenSize = [60,18];
    set(gca,'FontSize',font_size,'LineWidth',line_width_axis,'FontName','Arial',...
        'XTick',[0,T],'XTickLabel',{'$0$','$T$'},'YTick',[0,1],'TickLabelInterpreter',...
        font_interpreter,'Layer','Bottom','Color',1.0*ones(1,3))
    text(0.52*T,0.89,['Case ',Case],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = erre(d);
    text(0.52*T,0.78,['$\varepsilon_{\mathrm{e}} = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),...
        '\times 10^{',num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = errw(d);
    text(0.52*T,0.67,['$\varepsilon_{\mathrm{w}} = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),...
        '\times 10^{',num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = lambda;
    text(0.52*T,0.56,['$\lambda = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),...
        '\times 10^{',num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = mu;
    text(0.52*T,0.45,['$\mu = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),...
        '\times 10^{',num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = alpha;
    text(0.52*T,0.34,['$\alpha = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),...
        '\times 10^{',num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    val = T;
    text(0.52*T,0.23,['$T = ',num2str(val/10^(floor(log10(abs(val)))),'%.2f'),'\times 10^{',...
        num2str(floor(log10(abs(val))),'%i'),'}$'],'Interpreter','LaTeX',...
        'FontSize',font_size-2)
    box on
    drawnow
    if save_figs
        pause(1);
        print(gcf,[path_name,'Case',Case,num2str(d)],'-depsc2')
    end
    
    pause(0.1)
    
end
