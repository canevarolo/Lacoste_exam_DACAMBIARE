% Esame 2 febbraio 2024
% Simone Canevarolo
% S269893
% 27/06/2024

clear all
close all
clc

raggio = 20e-2; %m
alt = 7.5e3; % m
uu = 5; % m/s
T0 = 233; % K
hh = 100; % W/m*K
Taria = 270; % K
rovol = 1e3; % kg/m^3

cond = @(T) 5*(T/T0);
cp = @(T) 670+1.2*(T-233);

drvett = [5e-4, 1e-3, 5e-3, 1e-2]; % m
zz = 1;

dt = 1; % s
tfin = alt/uu;
tt = (0:dt:tfin);
Nt = length(tt);

for zz = 1:length(drvett)

    dr = drvett(zz);
    rr = (0:dr:raggio)';
    Nr = length(rr);

    for 

    Tm = T0*ones(Nr,1);
    Tcentro = T0*ones(Nt,1);
    Tsuperficie = T0*ones(Nt,1);

    aa = cond(Tm)*dt./rovol./cp(Tm)./dr^2;

    sup_diag = -aa.*ones(Nr,1);
    main_diag = (1+2*aa).*ones(Nr,1);
    sub_diag = -aa.*ones(Nr,1);
    
    Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

    AA = spdiags(Band,-1:1,Nr,Nr);

    bb = Tm;

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(end,end-1) = cond(Tm(end))/dr;
    AA(end,end) = -cond(Tm(end))/dr-hh;
    bb(end) = -Taria*hh;

    TT = AA\bb;

    err(zz) = norm(TT-Tm)/norm(Tm-T0);

    Tm = TT;

end

figure(1)
loglog(drvett,err,'LineWidth',2)
hold on