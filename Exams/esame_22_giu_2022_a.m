% Esame 15 giugno 2022
% Simone Canevarolo
% S269893
% 30 giugno 2024

clear all
close all
clc

raggio = 30e-2; % m
rovol = 917; % kg/m^3
cp = 1050; % J/kg/K
kk = 2.5; % W/m/K
Tini = 250; % K
press = 1e5; % Pa
hh = 25; % W/m^2/K

tfin = 3600; % s
Telio = @(t) 250+50/3600*t;

dr = 1e-3; % m
rr = (0:dr:raggio)';
Nr = length(rr);

dt = 1; % s
tt = (0:dt:tfin);
Nt = length(tt);

aa = kk*dt/dr^2/rovol/cp;

sub_diag = -aa*(1-dr./rr);
main_diag = (1+2*aa)*ones(Nr,1);
sup_diag = -aa*(1+dr./rr);

% sub_diag = aa*ones(Nr,1);
% main_diag = -2*aa*ones(Nr,1);
% sup_diag = aa*ones(Nr,1);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,Nr,Nr);
Tm = Tini*ones(Nr,1);
bb = Tm;

Tcentro = Tini*ones(Nt,1);
Tsuperficie = Tcentro;

for ii = 2:Nt

    Tout = Telio(ii);
    bb = Tm;

    AA(1,1) = 1;
    AA(1,2) = -1;
    bb(1) = 0;

    AA(end,end-1) = -kk/dr;
    AA(end,end) = kk/dr+hh;
    bb(end) = Tout*hh;

    TT = AA\bb;

    Tcentro(ii) = TT(1);
    Tsuperficie(ii) = TT(end);

    Tm = TT;

end

figure(1)
plot(tt,Tcentro)
hold on
plot(tt,Tsuperficie)
title('Andamento nel tempo della temperatura in centro e superficie')
xlabel('tempo [s]')
ylabel('Termperatura [K]')

figure(2)
plot(rr,TT)
title('Andamento del transitorio a t=3600 s')
xlabel('raggio [m]')
ylabel('Temperatura [K]')


Tplotzero = Tini*ones(Nt,1);

As = 4*pi*raggio^2;
Vol = 4/3*pi*raggio^3;

for ii = 2:Nt

    cc = hh*As/Vol*dt/rovol/cp;

    Tzerod = (-cc*Telio(ii)+Tm)/(1-cc);

    Tplotzero(ii) = Tzerod;

    Tm = Tzerod;

end

figure(3)
plot(tt,Tplotzero)



