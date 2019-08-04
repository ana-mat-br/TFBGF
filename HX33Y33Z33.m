% HOW TO CITE THIS CODE
%------------------------------------
% A. P. Fernandes, M. B. dos Santos, and G. Guimarães. 
% An analytical transfer function method to solve inverse heat conduction problems. 
% Applied Mathematical Modelling, 2015. ISSN 0307-904X. doi: http://dx.doi.org/10.1016/j. apm.2015.02.012.
%------------------------------------

tic
clc 
close all
clear all


prop = load('propX33Y33Z33.txt')';
kCOND=prop(1,:);
alfa=prop(2,:);
L=prop(3,:);
W=prop(4,:);
R=prop(5,:);

% time discretization
t=load('timeX33Y33Z33.txt')';  
dt=t(2)-t(1);  

% position
P=load('postionX33Y33Z33.txt')';  
x=P(1,:);
y=P(2,:);
z=P(3,:);

% heat flux area
areafluxo=load('areafluxoX33Y33Z33.txt')';
L1=areafluxo(1,:);
L2=areafluxo(2,:);
R1=areafluxo(3,:);
R2=areafluxo(4,:);

% convection
conv = load('convX33Y33Z33.txt')';
h1=conv(1,:);
h2=conv(2,:);
h3=conv(3,:);
h4=conv(4,:);
h5=conv(5,:);
h6=conv(6,:);
% Biot
B1=h1*L/kCOND;
B2=h2*L/kCOND;
B3=h3*W/kCOND;
B4=h4*W/kCOND;
B5=h5*R/kCOND;
B6=h6*R/kCOND;
% eigenvalues 
m=10;
n=10;
p=10;
alfam=abs(eigenvaluesX33(m,B1,B2));
betan=abs(eigenvaluesX33(n,B3,B4));
gamap=abs(eigenvaluesX33(p,B5,B6));

for s=2:length(t)

    somatorio=0;

    for k=1:p           
        for j=1:n
            for i=1:m           

                Amnp =  ( ( alfam(i)/L )^2 + ( betan(j)/W )^2 + ( gamap(k)/R )^2 ) * alfa;

                parcela = ...
                exp( - Amnp * t(s) ) ...
                * ( ( ( alfam(i) * cos( alfam(i) * x/L ) + B1 * sin( alfam(i) * x/L ) )  / ( ( alfam(i)^2 + B1^2 ) * ( 1 + ( B2 / ( alfam(i)^2 + B2^2 ) ) ) + B1 ) ) ) ...
                * ( ( ( betan(j) * cos( betan(j) * y/W ) + B3 * sin( betan(j) * y/W ) )  / ( ( betan(j)^2 + B3^2 ) * ( 1 + ( B4 / ( betan(j)^2 + B4^2 ) ) ) + B3 ) ) ) ...
                * ( ( betan(j) * cos( betan(j) ) + B3 * sin( betan(j) ) ) ) ...
                * ( ( ( gamap(k) * cos( gamap(k) * z/R ) + B5 * sin( gamap(k) * z/R  ) ) / ( ( gamap(k)^2 + B5^2 ) * ( 1 + ( B6 / ( gamap(k)^2 + B6^2 ) ) ) + B5 ) ) ) ...
                * (  sin( alfam(i) * L2/L )  -  sin( alfam(i) * L1/L )  -   ( B1 / alfam(i) ) * ( cos( alfam(i) * L2/L ) ) +  ( B1 / alfam(i) ) * ( cos( alfam(i) * L1/L ) ) )  ...
                * (  sin( gamap(k) * R2/R )  -  sin( gamap(k) * R1/R )  -   ( B5 / gamap(k) ) * ( cos( gamap(k) * R2/R ) ) +  ( B5 / gamap(k) ) * ( cos( gamap(k) * R1/R ) ) ) ;
                somatorio = somatorio + parcela; 

            end
        end
    end

    H(s) = (alfa/kCOND)*(8/W) *  somatorio;

end



jump=1;
figure(1)
plot(t(1:jump:length(H)),H(1:jump:length(H)),'-k')
grid
xlabel('time (s)')
ylabel('impulsive response [K/Js]')

save('../5/HX33Y33Z33.txt', 'H', '-ASCII')

