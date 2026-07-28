% Esame del 2 febbraio 2024 punto b - turno mattino 
% Simone Canevarolo
% S269893
% 19 aprile 2024

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

tempomax = altezza/uu; % s

dtvett = [0.5 1 5 10 50 100 500];

Tm = T0*ones(Nr,1);
zz = 1;

err = zeros(length(dtvett),1);

% Tcentro = T0*ones(Nt,1);
% Tsuperficie = T0*ones(Nt,1);


for zz = 1:length(dtvett)

    dt = dtvett(zz);
    tt = (0:dt:tempomax);
    Nt = length(tt);
    

for ii = 1:Nt

    aa = kk(Tm)*dt./dr^2./ro./cp(Tm);
    
%     sub_diag = -aa.*(1/dr^2-1/dr./rr);
%     main_diag = (1+2*aa);
%     sup_diag = -aa.*(1/dr^2-1/dr./rr);
    
    sub_diag = -aa;
    main_diag = (1+2*aa);
    sup_diag = -aa;
    
    Band = [sub_diag, main_diag, sup_diag];
    
    AA = spdiags(Band,-1:1,Nr,Nr);
    
    bb = Tm;

    AA(1,1) = -1;
    AA(1,2) = 1;
    bb(1) = 0;

    AA(end,end-1) = -kk(Tm(end))/dr;
    AA(end,end) = kk(Tm(end))/dr-hh;
    bb(end) = -Taria*hh;

    TT = AA\bb;

    err(zz) = norm(TT-Tm)/norm(Tm-T0);

    Tm = TT;

end

end

figure(1)
loglog(dtvett,err,'linewidth',2)
grid on
xlabel('Intervallo di discretizzazione [s]')
ylabel('Errore relativo')
title('Andamento dello errore relativo')

