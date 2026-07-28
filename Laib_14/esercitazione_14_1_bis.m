% Esercitazione 14, esercizio 1
% Simone Canevarolo
% S269893
% 24/02/2024

clear all
close all
clc

dd = 1e-2; % diametro interno, m
ss = 1e-3; % spessore, m
portata = 0.05; % portata, kg/s
ll = 2; % lunghezza, m
q0 = 3e3; % generazione di calore, W/m^3

Qv = @(x) q0*exp(-5.*(x-1).^2);

Tx0 = 25+273; % temperatura di ingresso, K
hh = 5e3; % coefficiente di scambio termico, W/m^2/K

roacc = 7.8e3; % densità, kg/m^3
cond = 30; % conducibilità termica, W/m/K
cpacc = 5e2; % calore specifico, J/kg/K

tfin = 10; % tempo, s
dt = 0.5;
tt = (0:dt:tfin);
Nt = length(tt);

ttvett = [1 3 5 10];
zz = 1;

dx = 1e-2; % m
xx = (0:dx:ll)';
Nx = length(xx);

Abase = pi*(dd/2)^2; % m^2
As = 2*pi*(dd/2)*ll; % m^2
Vol = Abase*ll; % m^3

uu = portata/Abase; % velocità fluido, m/s

Tm = Tx0*ones(2*Nx,1);
ii = 1;

for ii = 1:tfin

    % ACCIAIO

    aa1 = dt*cond/roacc/cpacc/dx^2;
    aa2 = hh*Abase*dt/Vol/roacc/cpacc;
    aa3 = dt/roacc/cpacc;

    % FLUIDO

    bb1 = uu*dt/dx;
    bb2 = hh*Abase*dt/roacc/cpacc/Vol;

    qvol = Qv(xx)/As;

    subsub_diag = zeros(2*Nx,1);
    sub_diag = zeros(2*Nx,1);
    main_diag = zeros(2*Nx,1);
    sup_diag = zeros(2*Nx,1);
    supsup_diag = zeros(2*Nx,1);

    % Fluido

    subsub_diag(1:2:end) = (-aa1)*ones(Nx,1);
    subsub_diag(2:2:end) = (-bb1)*ones(Nx,1);

    % sub_diag(1:2:end) = zeros(Nx,1);
    sub_diag(2:2:end) = (-bb2)*ones(Nx,1);

    main_diag(1:2:end) = (1-2*aa1-aa2)*ones(Nx,1);
    main_diag(2:2:end) = (1+bb1+bb2)*ones(Nx,1);

    sup_diag(1:2:end) = (-aa2)*ones(Nx,1);
    % sup_diag(2:2:end) = zeros(Nx,1);

    supsup_diag(1:2:end) = (-aa1)*ones(Nx,1);
    % supsup_diag(2:2:end) = zeros(Nx,1);

    Band = [subsub_diag, sub_diag, main_diag, sup_diag, supsup_diag];

    AA = spdiags(Band,-2:2,2*Nx,2*Nx);

    BB = zeros(2*Nx,1);
    BB(1:2:end,1) = aa3*qvol;

    bb = Tm+BB;

    % Acciaio
    % Adiabatica
    AA(1,1) = 1;
    AA(1,3) = -1;  

    AA(end-1,end-3) = -1;
    AA(end-1,end-1) = 1;

    bb(1) = 0;
    bb(end-1) = 0;

    % Fluido
    % Dirichlet
    AA(2,2) = 1;

    bb(2) = Tx0;

    TT = AA\bb;

    TTsolido = TT(1:2:end);
    TTfluido = TT(2:2:end);

    if ii == ttvett(zz)
    
        figure(1)
        plot(xx,TTsolido,'linewidth',2)
        hold on

        figure(2)
        plot(xx,TTfluido,'LineWidth',2)
        hold on
    
        zz = zz+1;
    
    end
    
Tm = TT;

end