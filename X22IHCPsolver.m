% HOW TO CITE THIS CODE
%------------------------------------
% A. P. Fernandes, M. B. dos Santos, and G. Guimarães. 
% An analytical transfer function method to solve inverse heat conduction problems. 
% Applied Mathematical Modelling, 2015. ISSN 0307-904X. doi: http://dx.doi.org/10.1016/j. apm.2015.02.012.
%------------------------------------

clc
close all
clear all

tic

H=load('H.txt');
T=load('T.txt');
q=load('q.txt');
t=load('time.txt');

T0=T(1,1);

% time discretization
dt=t(2)-t(1);


% Temperature position 1 and impulsive response 1
H1=H(:,1);
T1=T(:,1)-T0;

% Temperature position 2 and impulsive response 2
H2=H(:,2);
T2=T(:,2)-T0;

% Temperature position 3 and impulsive response 3
H3=H(:,3);
T3=T(:,3)-T0;



NR=2^20;

% Heat flux estimated by T1 e H1 

Hfreq1=fft(H1,NR);
Tfreq1=fft(T1,NR);
qfreq1=(Tfreq1./Hfreq1);
qtime1=(ifft(qfreq1)/(dt));

% Heat flux estimated by T2 e H2
H2=flipud(H2);
Hfreq2=fft(H2,NR);
Tfreq2=fft(T2,NR);
qfreq2=(Tfreq2./Hfreq2);
qtime2=(ifft(qfreq2)/(dt));

% Heat flux estimated by T3 e H3 
H3=flipud(H3);
Hfreq3=fft(H3,NR);
Tfreq3=fft(T3,NR);
qfreq3=(Tfreq3./Hfreq3);
qtime3=(ifft(qfreq3)/(dt));

%Table=[qtime1(1:length(t)) qtime2(1:length(t)) qtime3(1:length(t))];

jump=1;
ti=1;
tf=length(t)-50;

figure(1)
plot(t(ti:jump:tf),q(ti:jump:tf),'-k', 'LineWidth',2)
hold on
plot(t(ti:jump:tf),qtime1(ti:jump:tf), '-or', 'MarkerSize',8)
hold on
plot(t(ti:jump:tf),qtime2(ti:jump:tf), '-og', 'MarkerFaceColor','k', 'MarkerSize',3)
hold on
plot(t(ti:jump:tf),qtime3(ti:jump:tf), '-+b', 'MarkerSize',10)
xlabel('time (s)')
ylabel('heat flux (K/Js)')
legend ('q-sim', 'q1-est', 'q2-est','q3-est' ,'Location', 'Best')
set(gca, 'FontSize',14)
grid
hold off

toc