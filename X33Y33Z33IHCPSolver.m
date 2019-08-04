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


T=load('TX33Y33Z33.txt');
T0=25;
T=T-T0;
H=load('HX33Y33Z33.txt');

  
% time
t=load('timeX33Y33Z33.txt')';  
dt=t(2)-t(1); 
 
% heat flux simulated
q=load('qX33Y33Z33.txt')';

% IHCPSolver
NR=2^25;
Tfreq=fft(T,NR);
Hfreq=fft(H,NR);
qfreq=(Tfreq./Hfreq);
qtempo=ifft(qfreq)/dt;
qtempo=qtempo(1:length(q));


jump=1;
ti=1;
tf=length(t)-10;

figure(1)
plot(t(ti:jump:tf),q(ti:jump:tf),'-k','LineWidth',2)
hold on
plot(t(ti:jump:tf),qtempo(ti:jump:tf)','-or')
legend('q-sim', 'q-est', 'location', 'best')
xlabel('time (s)')
ylabel('heat flux (W/m^2)')
grid

