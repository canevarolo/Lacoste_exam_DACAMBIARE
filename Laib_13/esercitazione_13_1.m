% Esercitazione 13, esercizio 1
% Simone Canevarolo
% S269893
% 23/02/2024

clear all
close all
clc

ro = 140; % densità, kg/m^3
cp = 6e3; % calore specifico a pressione costante, J/kg/K
dd = 5e-3; % diametro, m
ll = 1; % lunghezza, m
uu = 2.5e-3; % velocità, m/s

dx = 1e-4;
xx = (0:dx:ll)';
Nx = length(xx);

Tin = 5; % temperatura iniziale del ferro, K
Tx0 = 10; % temperatura appena parte il tempo, K

tfin = 100; % tempo finale, s
dt = 1; % s
tt = (0:dt:tfin);
Ntt = length(tt);

tplot = [25 50 75 100];
iplot = round(tplot/dt)+1;
count = 1;

% Matrice

aa = uu*dt/dx;
main_diag = (1+aa)*ones(Nx,1);
sub_diag = -aa*ones(Nx,1);

Band = [sub_diag, main_diag];

AA = spdiags(Band,-1:0,Nx,Nx);

% Dirichlet all'inizio
AA(1,1) = 1;

TTin = Tin*ones(Nx,1);
Tguess = TTin;

figure(1)
plot(xx,TTin,'linewidth',3,'displayname','Temperatura iniziale')
legend
hold on

for ss = 1:Ntt
    
    Tguess(1) = Tx0;
    
    TT = AA\Tguess;
    Tguess = TT;

    if ss == iplot(count)

        plot(xx,TT,'LineWidth',2)
        count = count+1;

    end
    
    pause(0.1)

end

xlabel('Lunghezza')
ylabel('Temperatura')




