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

% thermophisical properties
kCOND=24;
alfa=7.0868e-06;
% geometry dimension x=L; y=W; z=R
L=1e-2;
W=1e-2;
R=10e-2;
prop = [kCOND alfa L W R];
save('../5/propX33Y33Z33.txt', 'prop', '-ASCII')

% initial temperature and temperature ambient
T0=25;
Tinf=30;
Teta0=T0-Tinf;

% time discretization
dt=1;
t=0:dt:200;
save('../5/timeX33Y33Z33.txt', 't', '-ASCII')

% calculated temperature, T(x,y,z) in the position 
x=3*L/5;
y=W;
z=3*R/50;
P=[x y z];
save('../5/postionX33Y33Z33.txt', 'P', '-ASCII')

% Simulated step heat flux 
c1=1e6;
c2=100;
q=c1*rectpuls(t-c2,c2);


% Simulated Gaussian heat flux 
% c1=1e7;
% c2=100;
% c3=1000;
% q= c1* (1/sqrt(2*pi))* exp (- ((t-c2).^2 )/c3 );

% constant heat flux
% q=1e8*ones(1,length(t));


save('../5/qX33Y33Z33.txt', 'q', '-ASCII')




% heat flux area
Li=0;   % initial value in x
Lf=L/5; % final value in x
Ri=0;   % initial value in z
Rf=R/50;% final value in z
L1=Li;
L2=Lf;
R1=Ri;
R2=Rf;
areafluxo=[L1 L2 R1 R2];
save('../5/areafluxoX33Y33Z33.txt', 'areafluxo', '-ASCII')

% convection coefficients
h=20;
h1=h;
h2=h;
h3=h;
h4=h; 
h5=h;
h6=h;
conv=[h1 h2 h3 h4 h5 h6];
save('../5/convX33Y33Z33.txt', 'conv', '-ASCII')

% Biot
B1=h1*L/kCOND;
B2=h2*L/kCOND;
B3=h3*W/kCOND;
B4=h4*W/kCOND;
B5=h5*R/kCOND;
B6=h6*R/kCOND;

% eigenvalues - m for x-direction; n for y-direction; p for z-direction 
m=10;
n=10;
p=10;
alfam=abs(eigenvaluesX33(m,B1,B2));
betan=abs(eigenvaluesX33(n,B3,B4));
gamap=abs(eigenvaluesX33(p,B5,B6));

for d=1:length(t)-1

    somatorioexp1=0;
    somatorioexp2=0;

    for k=1:p         
        for j=1:n
            for i=1:m
                
                % exponential term
                Amnp =  ( ( alfam(i)/L )^2 + ( betan(j)/W )^2 + ( gamap(k)/R )^2 ) * alfa;
                
                % temperature initial term
                parcelaexp1 = ...
                exp( - Amnp * t(d) ) ...
                * ( ( alfam(i) * cos( alfam(i) * x/L ) + B1 * sin( alfam(i) * x/L ) ) / ( ( ( alfam(i)^2 + B1^2 ) * ( 1 + ( B2 / ( alfam(i)^2 + B2^2 ) ) ) ) + B1 ) ) ...
                * ( ( betan(j) * cos( betan(j) * y/W ) + B3 * sin( betan(j) * y/W ) ) / ( ( ( betan(j)^2 + B3^2 ) * ( 1 + ( B4 / ( betan(j)^2 + B4^2 ) ) ) ) + B3 ) ) ...
                * ( ( gamap(k) * cos( gamap(k) * z/R ) + B5 * sin( gamap(k) * z/R ) ) / ( ( ( gamap(k)^2 + B5^2 ) * ( 1 + ( B6 / ( gamap(k)^2 + B6^2 ) ) ) ) + B5 ) ) ...
                * ( 1 / ( alfam(i) * betan(j) * gamap(k) ) ) ...
                * ( alfam(i) * sin( alfam(i) ) - B1 * ( cos( alfam(i) ) - 1  ) ) ...
                * ( betan(j) * sin( betan(j) ) - B3 * ( cos( betan(j) ) - 1  ) ) ...
                * ( gamap(k) * sin( gamap(k) ) - B5 * ( cos( gamap(k) ) - 1  ) );
                somatorioexp1 = somatorioexp1 + parcelaexp1;

                % heat flux term
                somatorioint = 0;

                for f=1:d
                    arg1 = ( t(d+1) - t(f+1) );
                    arg2 = ( t(d+1) - t(f) );
                    parcelaint = q(f) * (  exp(-Amnp*arg1) - exp(-Amnp*arg2)  );
                    somatorioint = somatorioint + parcelaint;
                end
                
                parcelaexp2 = ...
                  ( ( ( alfam(i) * cos( alfam(i) * x/L ) + B1 * sin( alfam(i) * x/L ) )  / ( ( alfam(i)^2 + B1^2 ) * ( 1 + ( B2 / ( alfam(i)^2 + B2^2 ) ) ) + B1 ) ) ) ...
                * ( ( ( betan(j) * cos( betan(j) * y/W ) + B3 * sin( betan(j) * y/W ) )  / ( ( betan(j)^2 + B3^2 ) * ( 1 + ( B4 / ( betan(j)^2 + B4^2 ) ) ) + B3 ) ) ) ...
                * ( ( betan(j) * cos( betan(j) ) + B3 * sin( betan(j) ) ) ) ...
                * ( ( ( gamap(k) * cos( gamap(k) * z/R ) + B5 * sin( gamap(k) * z/R  ) ) / ( ( gamap(k)^2 + B5^2 ) * ( 1 + ( B6 / ( gamap(k)^2 + B6^2 ) ) ) + B5 ) ) ) ...
                * (  sin( alfam(i) * L2/L )  -  sin( alfam(i) * L1/L )  -   ( B1 / alfam(i) ) * ( cos( alfam(i) * L2/L ) ) +  ( B1 / alfam(i) ) * ( cos( alfam(i) * L1/L ) ) )  ...
                * (  sin( gamap(k) * R2/R )  -  sin( gamap(k) * R1/R )  -   ( B5 / gamap(k) ) * ( cos( gamap(k) * R2/R ) ) +  ( B5 / gamap(k) ) * ( cos( gamap(k) * R1/R ) ) ) ...
                * (1/Amnp) * somatorioint;     
                
            
                somatorioexp2 = somatorioexp2 + parcelaexp2; 

            end
        end
    end

    T(d)=...
    (8*Teta0) ...
    *  somatorioexp1  ...
    +  (alfa/kCOND)*(8/W) ...
    *  somatorioexp2  + Tinf;

end

figure(1)
plot(t(1:length(q)),q(1:length(q)),'-k')
grid
xlabel('time (s)')
ylabel('heat flux (W/m^2)')

figure(2)
plot(t(1:length(T)),T(1:length(T)),'-k')
xlabel('time (s)')
ylabel('temperature (^oC)')
legend ('T(x,y,z)', 'Location', 'Best')
grid

save('../5/TX33Y33Z33.txt', 'T', '-ASCII')

toc