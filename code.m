clear all;
close all;
clc;

% Setup

L = 3;                  % Domain Length
n = 31;                 % No. of discrete points from source to sink
x = linspace(0,3,n);
deltax = x(2)-x(1);
gamma = 1.4;   

% Nozzle Shape
A = 1+2.2*(x-1.5).^2;   % Area
b = A.*(-1);            % For plotting the nozzle's shape's bottom portion
throat = find(A==1);
 
n_time = 1400;      % Number of time steps
C = 0.5;            % Courant number

% Initial Conditions
rho = 1-0.3146*x;  % Non dimensionnal Density rho/rho_0
T = 1-0.2314*x;    % Non dimensional Temperature T/T_0
v = (0.1+1.09*x).*T.^0.5;  % Non dimensional Velocity v/v_0 
p = rho.*T;  % Non dimensional Pressure P/P_0
    
time_nc = 0;

% Time Loop
for i = 1:n_time
    
    tic;
    dt_nc = min(C*(deltax./(v+(T.^(0.5)))));
    
    rho_old = rho;    
    v_old = v;   
    T_old = T;
    p_old = p;

    % Predictor method

    for j = 2:n-1

        dv_dx = (v(j+1)-v(j))/deltax;   
        drho_dx = (rho(j+1)-rho(j))/deltax;   
        dT_dx = (T(j+1)-T(j))/deltax;    
        dloga_dx = (log(A(j+1))-log(A(j)))/deltax;

        drho_dt_pred(j) = -rho(j)*dv_dx-rho(j)*v(j)*dloga_dx - v(j)*drho_dx;
        dv_dt_pred(j) = -v(j)*dv_dx-(1/gamma)*(dT_dx + (T(j)/rho(j))*drho_dx);
        dT_dt_pred(j) = -v(j)*dT_dx - ((gamma-1)*T(j))*(dv_dx + v(j)*dloga_dx);   

        % Temporary Values of parameters
        rho(j) = rho(j) + drho_dt_pred(j)*dt_nc;
        v(j) = v(j) + dv_dt_pred(j)*dt_nc;
        T(j) = T(j) + dT_dt_pred(j)*dt_nc;

    end


    % Corrector method  
    for j = 2:n-1  
        
        dv_dx = (v(j)-v(j-1))/deltax;
        drho_dx = (rho(j)-rho(j-1))/deltax; 
        dT_dx = (T(j)-T(j-1))/deltax;   
        dloga_dx = (log(A(j))-log(A(j-1)))/deltax;

        drho_dt_cor(j) = -rho(j)*dv_dx - rho(j)*v(j)*dloga_dx - v(j)*drho_dx; 
        dv_dt_cor(j) = -v(j)*dv_dx - (1/gamma)*(dT_dx + (T(j)/rho(j))*drho_dx);
        dT_dt_cor(j) = -v(j)*dT_dx - (gamma-1)*T(j)*(dv_dx + v(j)*dloga_dx);   

    end


    % Average time derivative  
    drho_dt_av = 0.5*(drho_dt_pred + drho_dt_cor);  
    dv_dt_av = 0.5*(dv_dt_pred + dv_dt_cor);   
    dT_dt_av = 0.5*(dT_dt_pred + dT_dt_cor);


    % Updating the flow field variables   
    for k = 2:n-1  
        rho(k) = rho_old(k) + drho_dt_av(k)*dt_nc;
        v(k) = v_old(k) + dv_dt_av(k)*dt_nc;      
        T(k) = T_old(k) + dT_dt_av(k)*dt_nc;  
    end


    % Boundary conditions Linear extrapolation using 2 nearest points    
    % Inlet    
    v(1) = 2*v(2) - v(3);    

    % Outlet    
    v(n) = 2*v(n-1) - v(n-2);    
    rho(n) = 2*rho(n-1) - rho(n-2);    
    T(n) = 2*T(n-1) - T(n-2);

    % Pressure    
    p = rho.*T;   
    
    % Finding the values of parameters at various values of timesteps
    if i==100
        flow_rate(1,:) = rho.*v.*A;
        flow_rate_nt(1) = i;
    elseif i==200
        flow_rate(2,:) = rho.*v.*A;
        flow_rate_nt(2) = i;
    elseif i==300
        flow_rate(3,:) = rho.*v.*A;
        flow_rate_nt(3) = i;
    elseif i==500
        flow_rate(4,:) = rho.*v.*A;
        flow_rate_nt(4) = i;
    elseif i==700
        flow_rate(5,:) = rho.*v.*A;
        flow_rate_nt(5) = i;
    elseif i==1400
        flow_rate(6,:) = rho.*v.*A;
        flow_rate_nt(6) = i;
    end    
end

    
% Plotting the result

figure();
    plot(x, A,'k', LineWidth=1.5);
    title('Nozzle Shape');
    hold on;
    plot(x, b, 'k', LineWidth=1.5);
    plot([x(round(n/2)) x(round(n/2))],[-5 5]);
    plot([x(1) x(n)],[0 0])

figure();
    plot(x,p);
    title('Pressure');
    grid on;
    xlabel('x / L');    ylabel('P/Po');

figure();
    plot(x,rho);
    title('Density');
    grid on;
    xlabel('x / L');    ylabel('ρ/ρo');
    
figure();
    plot(x,T);
    title('Temperature');
    grid on;
    xlabel('x / L');    ylabel('T/To');

figure();
    plot(x,v);
    title('Velocity');
    grid on;
    xlabel('x / L');    ylabel('V/V_0');

 
% Plotting Non-Dimensional Mass Flow rate at different No. of Time Steps
figure4 = figure('position',[50 -30 640 480]);
for q = 1:6
    plot(x,flow_rate(q,:),'linewidth',1.5);
    hold on;
end

% Legend at different time steps
title('Mass Flow rate at different No. of Time Steps');
grid on;    hold on;    xlabel('x / L');    ylabel('ρ*A*V / ρo*A*Vo');
n_steps1 = sprintf('%d Time Steps',flow_rate_nt(1));
n_steps2 = sprintf('%d Time Steps',flow_rate_nt(2));
n_Steps3 = sprintf('%d Time Steps',flow_rate_nt(3));
n_steps4 = sprintf('%d Time Steps',flow_rate_nt(4));
n_steps5 = sprintf('%d Time Steps',flow_rate_nt(5));
n_steps6 = sprintf('%d Time Steps',flow_rate_nt(6));
legend(n_steps1,n_steps2,n_Steps3,n_steps4,n_steps5,n_steps6);