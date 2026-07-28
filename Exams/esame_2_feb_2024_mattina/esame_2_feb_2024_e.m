% Esame 2 febbraio 2024
% Simone Canevarolo
% S269893
% 28/06/2024

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

dr = 1e-3; % m
rr = (0:dr:raggio)';
Nr = length(rr);

dt = 1; % s
tfin = alt/uu;
tt = (0:dt:tfin);
Nt = length(tt);

Asup = 4*pi*raggio^2;
Vol = 4/3*pi*raggio^3;
Lc = Asup/Vol;

Tm = T0*ones(Nr,1);

for ii = 2:Nt

    aa = hh*Lc*dt/rovol/cp(Tm);
    bb = -Taria*aa;

    sup_diag = (1-aa);
    main_diag = (1+2*aa);
    
    Band = [main_diag, sup_diag];

    AA = spdiags(Band,0:1,Nr,Nr);
    bb = Tm;

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(end,end-1) = cond(Tm(end))/dr;
    AA(end,end) = -cond(Tm(end))/dr-hh;
    bb(end) = -Taria*hh;

    Tzero = AA\bb;

    Tm = Tzero;

end

figure(1)
plot(tt,Tzero,'LineWidth',2)
hold on

Tm = T0*ones(Nr,1);

for ii = 2:Nt

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

    Tmedia = mean(TT);

    Tm = TT;

end

plot(tt,Tmedia,'LineWidth',2);



