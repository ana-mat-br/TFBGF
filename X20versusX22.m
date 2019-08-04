clc
close all
clear all

tic

m=100; % number of iterations = number of eigenvalues

% thermophisical properties
% pvc
k=0.15;
alfa=1.1708E-07;
% cooper
% k=401; 
% alfa=117E-06;

% geometry dimention 
%L=5E-3;
L=10E-3; % L inf 


% time discretization
dt=.1;
tf=180;
t=dt:dt:tf;

% constant heat flux
q0=1E3;

q=q0*ones(1,length(t));

%     %fluxo pulso retangular
%         c1=q0;
%         c2=30;
%         q=[c1*rectpuls(t-(1.1)*c2,c2)];

% temperature calculed at the point x, T(x,t) 
x=0E-3;


for c=1:length(t)
    
    for a=1:length(x)
        
        % X20 (pag. 188 + ierfc def. pag 501 - 2. ed. Beck)
        % ---
        
        % Temperature
        TX20(a,c) = (q0/k)*sqrt(4*alfa*t(c)) * ((1/sqrt(pi)) * exp(-(x(a)/sqrt(4*alfa*t(c)))^2)-((x(a)/sqrt(4*alfa*t(c))))*erfc((x(a)/sqrt(4*alfa*t(c)))));
        % Impuslive response
        HX20(a,c) = (alfa/k)*(1/(sqrt(4*pi*alfa*t(c))))*2*exp(-((x(a)^2)/4*alfa*t(c)));
        
        % X22 
        % ---
        somaTX22=0;
        somaHX22=0;
        for j=1:m
            
            Bn=(j*pi/L)^2 * alfa ;
            parcelaTX22 =  cos(j*pi*x(a)/L) * (exp(-Bn*t(c)))/j^2; 
            somaTX22 = somaTX22 + parcelaTX22;
            parcelaHX22 =  cos(j*pi*x(a)/L) * (exp(-Bn*t(c)));
            somaHX22 = somaHX22 + parcelaHX22;
             
        end

        % Temperature
        TX22(a,c)  = (q0*L/k)*( (alfa*t(c)/L^2) + ( (1/2)*(x(a)/L)^2 )  - (x(a)/L) + (1/3) - ( (2/pi^2)*somaTX22 ) );
        % Impuslive response
        HX22(a,c)  =  alfa/(k*L)+(alfa*2)/(k*L)*somaHX22;
             
        
    end
end


TX22conv=conv(q,HX22)*dt;

TX20conv=conv(q,HX20)*dt;


% ------------------------------------------------------------------------
% gráfico: temperatura
% ------------------------------------------------------------------------
figure(1)
jump=100;
plot(t(1:jump:length(TX20)),TX20(1:jump:length(TX20)),'ok')
hold on
plot(t(1:jump:length(TX20)),TX22(1:jump:length(TX22)),'+r')
hold on
plot(t(1:jump:length(TX20)),TX22conv(1:jump:length(TX22)),'o--k')
hold on
plot(t(1:jump:length(TX20)),TX20conv(1:jump:length(TX22)),'+--r')
xlabel('time (s)')
ylabel('Temperature (ºC)')
legend ('TX20', 'TX22' , 'TX22conv', 'TX20conv','Location', 'Best')
title(strcat('temperature calculed at the point x=',num2str(x),', dt=',num2str(dt)))
set(gca, 'FontSize',14)
grid
hold off

%Table 

format long
[t(1:length(TX20))' TX20(1:length(TX20))' TX22(1:length(TX20))' TX22conv(1:length(TX20))'  TX20(1:length(TX20))'-TX22(1:length(TX20))'  TX22conv(1:length(TX20))'-TX22(1:length(TX20))']