% Esercitazione 13, esercizio 1
% Simone Canevarolo
% S269893
% 15/03/2024

clear all
close all
clc

ro = 140; % densità, kg/m^3
cp = 6e3; % calore specifico, J/kg/K
dd = 5e-3; % diamentro, m
ll = 1; % lunghezza, m

uu = 0.25e-2; % velocità elio, m/s
Tini = 5; % temperatura iniziale barra, K
T0 = 10; % temperatura all'inizio della barra a t>0, K

dx = 1e-4;
xx = (0:dx:ll)';
Nx = length(xx);

dt = 1; % tempo, s
aa = uu*dt/dx;

tfin = 100; % tempo finale, s
tt = (0:dt:tfin);
Nt = length(tt);

tplot = [25 50 75 100];
iplot = round(tplot/dt)+1;
count = 1;

TTin = Tini*ones(Nx,1);
Tguess = TTin;

% non serve sup_diag
main_diag = (1+aa)*ones(Nx,1);
sub_diag = -aa*ones(Nx,1);

Band = [sub_diag, main_diag];

AA = spdiags(Band,-1:0,Nx,Nx);

AA(1,1) = 1;

figure(1)
plot(xx,TTin,'LineWidth',2)
hold on

for iterazioni = 1:Nt

    Tguess(1) = T0;

    TT = AA\Tguess;
    Tguess = TT;

    if iterazioni == iplot(count)

        plot(xx,TT,'LineWidth',2);
        count = count+1;

    end

end