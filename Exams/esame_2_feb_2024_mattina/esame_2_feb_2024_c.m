% Esame del 2 febbraio 2024 punto c- turno mattino
% Simone Canevarolo
% S269893
% 11 aprile 2024

clear all
close all
clc

raggio = 20e-2; % m
altezza = 7.5e3; % m
uu = 5; % m/s
T0 = 233; % K
hh = 100; % W/m^2/K
Taria = 270; % K
ro = 1e3; % kg/m^3

kk = @(T) 5*(T/T0);
cp = @(T) 670+1.2*(T-233);

dr = 1e-3;
rr = (0:dr:raggio)';
Nr = length(rr);

dt = 1; % s
tempomax = altezza/uu; % s
tt = (0:dt:tempomax);

ttvett = [100 500 1000];

Tm = T0*ones(Nr,1);
zz = 1;

for ii = 2:length(tt)

    aa = kk(Tm)*dt./dr^2./ro./cp(Tm);
    
%     sub_diag = -aa.*(1/dr^2-1./dr./rr);
%     main_diag = (1+2*aa);
%     sup_diag = -aa.*(1/dr^2-1./dr./rr);
    
    sub_diag = -aa;
    main_diag = (1+2*aa);
    sup_diag = -aa;
    
    Band = [[sub_diag(2:end);0], main_diag, [0;sup_diag(1:end-1)]];
    
    AA = spdiags(Band,-1:1,Nr,Nr);
    
    bb = Tm;

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(end,end-1) = kk(Tm(end))/dr;
    AA(end,end) = -kk(Tm(end))/dr-hh;
    bb(end) = -Taria*hh;

    TT = AA\bb;

if ii == ttvett(zz)

    figure(1)
    plot(rr,TT,'LineWidth',2)
    hold on
    grid on

    zz = zz+1;

end

Tm = TT;

end