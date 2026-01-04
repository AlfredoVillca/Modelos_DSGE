%------------------------------------------------------------------------
% Modelo de solow estocastico con simulaciones de ciclos
%
% Alfredo Villca
% alfredovillca569@gmail.com
%-------------------------------------------------------------------------

clc; clear all; close all;

%Parametrizacion
A = 1;
n = 0.02;
delta = 0.1;
alpha = 0.36;
s = 0.2;
% Estado estacionario
k_ss = ((s*A)/(n+delta))^(1/(1-alpha));
y_ss = A*(k_ss^alpha);
c_ss = (1-s)*y_ss;
s_ss = s*y_ss;

sigma = 0.01; %desviacion estandar
Nobs =100;   % número de periodos a simular
Nsim = 2000;   % numero de simulaciones
epsilon = sigma*randn(Nobs, Nsim);

tiempo=1:1:Nobs;
%Creando vectores para almzacenar resultados

K=zeros(Nobs, Nsim);
S=zeros(Nobs, Nsim);
Y=zeros(Nobs, Nsim);
C=zeros(Nobs, Nsim);


K(1,1:Nsim) = (s*A/(n+delta))^(1/(1-alpha));  
Y(1,1:Nsim) = K(1,1)^alpha;
S(1,1:Nsim) = s*Y(1,1);
C(1,1:Nsim) =(1-s)*Y(1,1); 
%MK(1,1:Nsim) = (s*A/(n+delta))^(1/(1-alpha)); 

% Iteraciones simuladas
for i = 2:Nobs
    for j = 1:Nsim
        K(i,j) = s*A*exp(epsilon(i, j))*(K(i-1, j)^alpha)+ (1-delta - n)*K(i-1, j);
       % K(i,j) = ( (s*A*exp(epsilon(i, j)))/(1+n) )*(K(i-1, j)^alpha)+ ((1-delta)/(1+n))*K(i-1, j);
        Y(i,j) = (K(i-1, j)^alpha);
        S(i,j) = s*Y(i,j);
        C(i,j) = (1-s)*Y(i,j);
    end
end

%prueba de las funciones de K de nuestro modelo vs MacCandless
%figure
%plot(K(:,1), 'Color','black', 'LineWidth',1.5), hold on, 
%plot(MK(:,1), 'Color','blue', 'LineWidth',1.5), hold off,
%legend('Nuestro', 'MacCandless')
%title('Capital')


lwdt = 1.5;
%series de los errores
figure
plot(tiempo, epsilon(:,1)), hold on,
plot(tiempo, epsilon(:,2)),
plot(tiempo, epsilon(:,3))
legend('\epsilon_{t,1}', '\epsilon_{t,2}', '\epsilon_{t,3}')

% Series temporales artificiales MODELO NO LINEAL
figure
subplot(221), plot(tiempo, K(:,1), 'b','LineWidth',lwdt), hold on, 
subplot(222), plot(tiempo, Y(:,1), 'b', 'LineWidth',lwdt), hold on, 
subplot(223), plot(tiempo, C(:,1),'b', 'LineWidth',lwdt), hold on, 
subplot(224), plot(tiempo, S(:,1), 'b', 'LineWidth',lwdt), hold on, 
title('Ahorro')




figure
subplot(221), plot(tiempo, K(:,1), 'b','LineWidth',lwdt), hold on, 
    plot(tiempo, K(:,2), 'r', 'LineStyle',':','LineWidth',lwdt), 
    plot(tiempo, K(:,3),'g','LineWidth',lwdt),    title('Capital')
    legend({'k_{1t}','k_{2t}','k_{3t}', 'k_{4t}','k_{5t}'},'NumColumns',3),
subplot(222), plot(tiempo, Y(:,1), 'b', 'LineWidth',lwdt), hold on, 
    plot(tiempo, Y(:,2), 'r', 'LineStyle',':','LineWidth',lwdt), 
    plot(tiempo, Y(:,3), 'g','LineWidth',lwdt),    title('Producto') 
    legend({'y_{1t}','y_{2t}','y_{3t}', 'y_{4t}','y_{5t}'}, 'NumColumns',3), 
subplot(223), plot(tiempo, C(:,1),'b', 'LineWidth',lwdt), hold on, 
    plot(tiempo, C(:,2), 'r', 'LineStyle',':','LineWidth',lwdt),  
    plot(tiempo, C(:,3), 'g','LineWidth',lwdt),  title('Consumo')
    legend({'c_{1t}','c_{2t}','c_{3t}', 'c_{4t}','c_{5t}'}, 'NumColumns',3), 
subplot(224), plot(tiempo, S(:,1), 'b', 'LineWidth',lwdt), hold on, 
    plot(tiempo, S(:,2), 'r', 'LineStyle',':','LineWidth',lwdt), 
    plot(tiempo, S(:,3), 'g', 'LineWidth',lwdt),
    legend({'s_{1t}','s_{2t}','s_{3t}', 's_{4t}','s_{5t}'}, 'NumColumns',3), 
    title('Ahorro')
    
    
   
%Analisis Boostrapp


% Momentos estadísticos
eps_men = mean(epsilon(:,:));

figure
histogram(eps_men, 'facecolor',[0 0 1])




mean_K = mean(K(:,:));
mean_Y = mean(Y(:,:));
mean_C = mean(C(:,:));
mean_S = mean(S(:,:));


figure
subplot(221), histogram(mean_K, 'facecolor',[0 0 1]), hold on
xline(k_ss,'r','linestyle','--','LineWidth',2)
title('Capital')

subplot(222), histogram(mean_Y, 'facecolor',[0 0 1]), hold on
xline(y_ss,'r','linestyle','--','LineWidth',2)
title('Producto')

subplot(223), histogram(mean_C, 'facecolor',[0 0 1]), hold on
xline(c_ss,'r','linestyle','--','LineWidth',2)
title('Consumo')

subplot(224), histogram(mean_S, 'facecolor',[0 0 1]), hold on
xline(s_ss,'r','linestyle','--','LineWidth',2)
title('Ahorro')
