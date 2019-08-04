
close all
clear all

% pvc
k=0.152;
alfa=1.24E-07;

% cooper
% k=401;        
% alfa=117E-06;

% Inicial temperature T0
T0=0;

% Time discretization
dt=1;
tf=100;
t=[dt:dt:tf];

%Heat flux W/m^2.
q0=100;
q=q0*ones(1,length(t));

% Temperatura position in (0,inf)
x=[0 .001 .01]; % [T1 T2 T3]


% ------------------------------------------------------------------------
for c=1:length(t)

    % X20 (pag. 188 + ierfc def. pag 501 - 2. ed. Beck)
    for a=1:length(x)        
        TX20(a,c) = (q0/k)*sqrt(4*alfa*t(c)) * ((1/sqrt(pi)) * exp(-(x(a)/sqrt(4*alfa*t(c)))^2)-((x(a)/sqrt(4*alfa*t(c))))*erfc((x(a)/sqrt(4*alfa*t(c)))));
    end
   
end

T1=TX20(1,:);
T2=TX20(2,:);
T3=TX20(3,:);


figure(1)
plot(t(1:length(T1)),T1(1:length(T1)),'-r')
hold on
plot(t(1:length(T2)),T2(1:length(T2)),'-g')
hold on
plot(t(1:length(T3)),T3(1:length(T3)),'-b')
xlabel('time (s)')
ylabel('Temperature (ºC)')
legend ('T1', 'T2', 'T3' ,'Location', 'Best')
set(gca, 'FontSize',14)
grid
hold off

