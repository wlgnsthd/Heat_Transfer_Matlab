%% Aerospace Engineering 20171648 Song Jihoon

%{
2D Steady State Heat Conduction ( i = space, k = time ) 
Assumption : Fixed boundaries temperature, Interior node, No heat
generation, n = m, T_0(m,n) = 100 (Celsius)

From heat diffusion equation,
dT = alpha*( T_(m+1,n)+T_(m-1,n)+T_(m,n+1)+T_(m,n-1)-4*T_(m,n) )*dt/(dx)^2

%}

%% Initialization
clear 
close all 
clc
tic % Start CPU Run Time Check 

%% Declaration of Constant Variable  
n = 52; % number of x index (divided into n-1 parts) max 52
L = 0.3; % Length of the bar[m]
alp = 8.5 * 10^(-5); % diffusivity
dx = L / (n-1); % difference of lenght[m]
dt = 0.1; % difference of time[s]

T_k = zeros(n,n) + 100; 
T_kk = zeros(n,n) + 100; 
T_k(:,1) = 0; 
T_kk(:,1) = 0;
T_k(:,n) = 0;
T_kk(:,n) = 0;
T_k(1,:) = 200;
T_kk(1,:) = 200;
T_k(n,:) = 200;
T_kk(n,:) = 200;
dT = zeros(n:n);
%% Initial Calculator
for i = 2:n-1 % except base tip (constant temp)
    for k = 2:n-1
        dT(i,k) = alp * dt / (dx)^2 * (T_k(i+1,k)+T_k(i-1,k)+T_k(i,k+1)+T_k(i,k-1)-4*T_k(i,k)); 
    end
end
T_kk = T_k + dT;
%% Iteration
avg_relative_err = sum(abs(dT),'all') / 100; % except base tip
passed_dt = 2;
while avg_relative_err >= 10^(-5)
    
    T_k = T_kk;
   
    for i = 2:n-1 % except base tip (constant temp)
        for k = 2:n-1
            dT(i,k) = alp * dt / (dx)^2 * (T_k(i+1,k)+T_k(i-1,k)+T_k(i,k+1)+T_k(i,k-1)-4*T_k(i,k)); 
        end
    end
   
    avg_relative_err = sum(abs(dT),'all') / 100; % except base tip

    T_kk = T_k + dT;

    passed_dt = passed_dt + 1;

end


%% Plot 1 (Question 1)
figure(1)
x = linspace(0,30,n);
y = linspace(0,30,n);
contourf(x,y,T_kk,'--','ShowText','on')
colormap('jet')
colorbar
title 'T(x,y) when the plate reaches steady state';
xlabel 'x [cm]';
ylabel 'y [cm]';


%% Plot 2 (Question 2)
figure(2)
x = linspace(0,30,n);
y = T_kk(n/2,:);
plot(x,y);
title 'T(x) @ y = 30 [cm] when the plate reaches steady state';
xlabel 'x [cm]';
ylabel 'T [Celsius]';

%% Elapsed Real Time
x = ['Elapsed Real Time = ',num2str(passed_dt*dt),' second(s)'];
disp(x)

%% Elapsed Real Time
x = ['CPU Run Time = ', num2str(toc),' second(s)'];
disp(x)