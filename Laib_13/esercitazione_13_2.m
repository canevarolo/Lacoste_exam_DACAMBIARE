% Esercitazione 13, esercizio 2
% Simone Canevarolo
% S269893
% 24/02/2024

clear all
close all
clc

dd = 5e-3; % diametro tubo, m
ll = 1; % lunghezza tubo, m
Tin = 5; % temperatura iniziale tubo, K
Tx0 = 10; % temperatura all'ingresso fluido, K

tfin = 1; % s

% T = linspace(0,30,50);
ro = @(T) 10./(1+T);
% plot(T,ro(T))

dx = 1e-3;
xx = (0:dx:ll)';
Nx = length(xx);

dt = 1e-2; % s
tt = (0:dt:tfin);
Ntt = length(tt);
TTin = Tin*ones(Nx,1);

% aa = ro(Tin)*dt/dx*ones(Nx,1);
% count = 1;
% 
% main_diag = 1+aa;
% sub_diag = -aa;
% 
% Band = [sub_diag, main_diag];
% AA = spdiags(Band,-1:0,Nx,Nx);
% 
% % condiziioni al contorno
% AA(1,1) = 1;

figure(1)
plot(xx,TTin);
hold on

Tm = TTin;
tplot = [0.25 0.5 0.75 1];
iplot = round(tplot/dt)+1;
% Nt = length(tplot);
count = 1;

for ii = 2:Ntt

    aa = ro(Tm)*dt/dx;
    
    main_diag = 1+aa;
    sub_diag = -aa;
    
    Band = [[sub_diag(2:end);0], main_diag];
    AA = spdiags(Band,-1:0,Nx,Nx);
    
    % condizioni al contorno
    AA(1,1) = 1;
    Tm(1) = Tx0;

    TT = AA\Tm;
    Tm = TT;

    if ii == iplot(count)

    plot(xx,TT)
    count = count+1;
    
    end

end


