% Esercitazione 13, esercizio 2
% Simone Canevarolo
% S269893
% 28/03/2024

clear all
close all
clc

dd = 5e-3; % diametro, m
ll = 1; % lunghezza, m
% ll >> dd, applico ipotesi 1D

ro = @(T) 10./(1+T);
% T = (0:100);

% figure(1)
% plot(T,ro(T))

Tiniz = 5; % temperatura all'istante iniziale lungo tutta la barra, K
T0 = 10; % temperatura della barra a x=0 lungo tutto il transitorio, K

dx = 1e-3;
xx = (0:dx:ll)';
Nx = length(xx);

% aa = ro*dt/dx;

tmax = 1;
dt = 0.01; % s
tt = (0:dt:tmax);
Nt = length(tt);

tplot = [0.10 0.20 0.50 0.75 1];
iplot = round(tplot/dt)+1;
count = 1;

Tguess = Tiniz*ones(Nx,1);

figure(1)
plot(xx,Tguess,'LineWidth',2)
hold on


for ii = 2:Nt

    aa = ro(Tguess)*dt/dx;

    main_diag = (1+aa);
    sub_diag = -aa;

    Band = [[sub_diag(2:end);0],main_diag];

    AA = spdiags(Band,-1:0,Nx,Nx);

    AA(1,1) = 1;
    Tguess(1) = T0;

    TT = AA\Tguess;
    Tguess = TT;

    if ii == iplot(count)

        plot(xx,TT,'LineWidth',2);
        count = count+1;

    end

end
