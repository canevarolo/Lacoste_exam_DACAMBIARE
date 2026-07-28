% Esercitazione 13, esercizio 3
% Simone Canevarolo
% S269893
% 24/02/2024

clear all
close all
clc

ro = 140; % densità, kg/m^3
cp = 6e3; % calore specifico a pressione costante, J/kg/K
dd = 5e-3; % diametro, m
ll = 1; % lunghezza, m
uu = 1e-2; % velocità, m/s
tfin = 100; % tempo finale, s

dx = 1e-4;
xx = (0:dx:ll)';
Nx = length(xx);

dt = 0.5; % s
tt = (0:dt:tfin);
Nt = length(tt);

aa = uu*dt/dx;

main_diag = (1+aa)*ones(Nx,1);
sub_diag = -aa*ones(Nx,1);

Band = [sub_diag, main_diag];

AA = spdiags(Band,-1:0,Nx,Nx);

% Condizioni al contorno, Robin in ingresso

bb = 
BB = 

AA(1,1) = 1;

tplot = [5 10 25 50 100];
iplot = round(tplot/dt)+1;
count = 1;


for ii = 2:tfin



    TT = AA\Tm;

    Tm = TT;


end




