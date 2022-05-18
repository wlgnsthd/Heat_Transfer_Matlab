%% Aerospace Engineering 20171648 Song Jihoon

%{
1D Unsteady State Heat Conduction ( i = space, k = time ) 
Assumption : 1D Heat Transfer, Constant Base Tip Temperature 

Laplace(T_i) = d^2 T / dx^2= (T_(i+1) -2 T_(i) + T_(i-1)) / (dx)^2
T_(k+1) = (alpha)(Laplace(T_i))(dt) + T_(k)
%}

%% Initialization
clear 
close all 
clc
tic % Start CPU Run Time Check 

%% Declaration of Constant Variable  
n = 50; % number of x index (divided into n-1 parts)
L = 0.3; % Length of the bar[m]
alp = 8.5 * 10^(-5); % diffusivity
dx = L / (n-1); % difference of lenght[m]
dt = 0.1; % difference of time[s]

T_k = zeros(1,n) + 500; % T(k,i)
T_kk = zeros(1,n); % T(k+1,i)
T_k(1,1) = 300; % base
T_kk(1,1) = 300;
T_k(1,n) = 400; % tip
T_kk(1,n) = 400;
T_result = T_k; % T(k=1,x)

lap_op = zeros(1,n); % laplace operator array

%% Initial Calculator
for i = 2:n-1 % except base tip (constant temp)
    lap_op(1,i) = (T_k(1,i+1) - 2*T_k(1,i) + T_k(1,i-1)) / (dx)^2;
end
 
for i = 2:n-1 % except base tip 
   T_kk(1,i) = dt *  alp * lap_op(1,i) + T_k(1,i);
end

%% Iteration
diff_avg = sum(abs(T_kk - T_k)) / (n - 2); % except base tip

k = 2; % start with k = 2

while diff_avg >= 10^(-6)
    for c = 1 : n
    T_result(k,c) = T_kk(1,c); % T(t,x)
    end

    T_k = T_kk;
   
    for i = 2:n-1 % except base tip
        lap_op(1,i) = (T_k(1,i+1) - 2*T_k(1,i) + T_k(1,i-1)) / (dx)^2;
    end
    
    for i = 2:n-1 % except base tip 
        T_kk(1,i) = dt *  alp * lap_op(1,i) + T_k(1,i);
    end
   
    diff_avg = sum(abs(T_kk - T_k)) / (n - 2); % except base tip

    k = k + 1;
end

%% Plot 1 (Question 1)
figure(1)
x = linspace(0,30,n);
y1 = T_result(1801,:); % 3min=180sec=(k=1801)
y2 = T_result(4801,:); % 8min=480sec=(k=4801)
plot(x,y1,x,y2);
title 'T(x) @ t = 3 & 8 [min]';
xlabel 'x [cm]';
ylabel 'T [K]';
legend('3min','8min');

%% Plot 2 (Question 2)
figure(2)
t = linspace(0,(k-1)*dt,k-1);
y = T_result(:,n/2);
plot(t,y);
title 'T(t) @ x = 15 [cm]';
xlabel 't [sec]';
ylabel 'T [K]';

%% Plot 3 (Question 3)
figure(3)
tiledlayout(2,1)
nexttile
x = linspace(0,30,n);
y = T_result(k-1,:);
plot(x,y);
title 'T(x) when the rod reaches steady state';
xlabel 'x [cm]';
ylabel 'T [K]';

nexttile
imagesc(y)
axis off
colorbar

%% Question 4
x = ['Elapsed Real Time = ',num2str(k*dt),' second(s)'];
disp(x)

%% Question 5
x = ['CPU Run Time = ', num2str(toc),' second(s)'];
disp(x)