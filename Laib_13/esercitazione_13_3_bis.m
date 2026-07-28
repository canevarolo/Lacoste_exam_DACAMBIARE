% Esercitazione 13, esercizio 3
% Simone Canevarolo
% S269893
% 30/03/2024

clear all
close all
clc

ro = 140; % densità, kg/m^3
cp = 6e3; % calore specifico, J/kg/K
dd = 5e-3; % diamentro, m
lltot = 1; % lunghezza, m
ll = lltot/2;

uu = 0.1e-2; % velocità elio, m/s
Tini = 5; % temperatura iniziale barra, K
T0 = 10; % temperatura all'inizio della barra a t>0, K

The = 5; % temperatura elio refrigerante, K
hhhe = 100; % coefficiente di scambio termico, W/m^2/K

dx = 1e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

tfin = 100; % tempo finale, s
dt = 0.1; % s
tt = (0:dt:tfin);
Nt = length(tt);

tplot = [5 10 25 50 100];
iplot = round(tplot/dt)+1;
count = 1;

Sup_base = pi*(dd/2)^2;
As = 2*pi*(dd/2)*ll;
Vol = Sup_base*ll;


Tm = Tini*ones(Nx,1);





