%-------------------------------------------------------------------------
% Simulation of the stochastic Solow model
%
%          Alfredo Villca
%     alfredovillca569@gmail.com
%--------------------------------------------------------------------------

clear all; close all; clc

% Parameters obtained from MacCandles
A         = 1;
n         = 0.02;
delta     = 0.1;
alpha     = 0.36;
s         = 0.2;
rho       = 0.9;
sigma_ee  = 0.01;

% Steady state
k_ss = ((s*A)/(n+delta))^(1/(1-alpha));
y_ss = A*(k_ss^alpha);
c_ss = (1-s)*y_ss;
s_ss = s*y_ss;

lwdt   = 2;
T      = 100;
tiempo = [1:1:T];

k  = zeros(1,T);
y  = zeros(1,T);
c  = zeros(1,T);
ss = zeros(1,T);

k(1)  = k_ss;
y(1)  = y_ss;
c(1)  = c_ss;
ss(1) = s_ss;

e(1) = 0;  % valor inicial del proceso estocástico

for t=2:T
    if t==2                     % periodo en la que se aplica el choque unico
        u = sigma_ee;           % magnitud del choque
        e(t) = rho*e(t-1) + u;  
    else                     
        u = 0;
        e(t) = rho*e(t-1) + u;
    end
k(t) = ((s*A)/(1+n))*(exp(e(t)))*(k(t-1)^alpha) + ((1-delta)/(1+n))*k(t-1);
y(t) = k(t-1)^alpha;
c(t) = (1-s)*y(t);
ss(t)= s*y(t);    
end

figure
plot(e, 'black', 'LineWidth',lwdt)
title('PTF')

figure
subplot(2,2,1), plot(tiempo, k,'black', 'LineWidth',lwdt), hold on, 
    plot(xlim,[k_ss k_ss],'black','LineStyle','--')
    xlabel('Periodo temporal'), ylabel('desviación del ss'), 
    title('Capital')

subplot(2,2,2), plot(tiempo, y, 'Color','black', 'LineWidth',lwdt), hold on,
    plot(xlim,[y_ss y_ss], 'color','black','LineStyle','--')
    xlabel('Periodo temporal'), ylabel('desviación del ss'), 
    title('Producto')

subplot(2,2,3), plot(tiempo, c, 'Color','black', 'LineWidth',lwdt), hold on,
    plot(xlim,[c_ss c_ss], 'color','black','LineStyle','--') 
    xlabel('Periodo temporal'), ylabel('desviación del ss'), 
    title('consumo')

subplot(2,2,4), plot(tiempo, ss, 'Color','black', 'LineWidth',lwdt), hold on,
    plot(xlim,[s_ss s_ss], 'color','black','LineStyle','--')     
    xlabel('Periodo temporal'), ylabel('desviación del ss'), 
    title('Ahorro')

