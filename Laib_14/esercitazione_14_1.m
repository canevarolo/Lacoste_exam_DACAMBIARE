% Esercitazione 14, esercizio 1
% Simone Canevarolo
% S269893
% 30/01/2024

clear all
close all
clc

%%

dd = 1e-2; % diametro interno, m
rin = dd/2; % m
spess = 1e-3; % spessore, m
ll = 2; % lunghezza, m
qzero = 3e3; % potenza iniziale, W/m
Tingr = 25; % temperatura ingresso acqua, °C
gg = 0.05 ; % portata, kg/s

roacc = 7800; % densità, kg/m^3
cpacc = 500; % calore specifico acciaio, J/Kg/K
kk = 30; % condicibilità, W/m/K

Q = @(x) qzero*exp(-5.*(x-1).^2);

tt = [1,3,5,10]; % tempo, s
hin = 5e3; % coefficiente di scambio termico, W/m^2/K

% Distribuzione a t = 10 s è stazionaria?

dx = 1e-2;
xx = (0:dx:ll)';
Nx = length(xx);

Ain = 2*pi*rin*ll; % m^2
Asez = pi*(rin+spess)^2-pi*(rin)^2; % Area sezione, m^2
Vf = pi*rin^2*ll; % m^3
Vs = pi*(rin+spess)^2*ll-Vf; % m^3

dt = 1e-1; % s

uu = gg/1000/Ain; % velocità nel condotto, m/s

aa = kk*dt/roacc/cpacc/dx^2;
bin = hin*Ain*dt/Vs/roacc/cpacc;
qq = Q(xx)*dt/roacc/cpacc/Asez;
cc = dt*uu/dx;
ee = hin*Ain*dt/1000/4186/Vf;

subsub = zeros(2*Nx,1);
sub = zeros(2*Nx,1);
main = zeros(2*Nx,1);
sup = zeros(2*Nx,1);
supsup = zeros(2*Nx,1);

subsub(1:2:end) = -aa;
subsub(2:2:end) = cc;

sub(1:2:end) = -ee;
sub(2:2:end) = 0;

main(1:2:end) = 1+2*aa+bin;
main(2:2:end) = 1+cc;

sup(1:2:end) = 0;
sup(2:2:end) = -bin;

supsup(1:2:end) = -aa;
supsup(2:2:end) = ee;

Band = [subsub,sub,main,sup,supsup];

AA = spdiags(Band,-2:2,2*Nx,2*Nx);

bb = zeros(2*Nx,1);
bb(1:2:2*Nx) = qq;
t1=1;
t2=3;
t3=5;
t4 = 10;

Tm = Tingr;
time=0;
iplot=0;

for counter = 1:t4
    time=time+dt;
    AA([1,2,end-1],:) = 0;
    BB = bb+Tm;
    AA(1,1) = 1;
    AA(1,3) = -1;
    BB(1) = 0;

    AA(2,2) = 1;
    BB(2) = Tingr;

    AA(end-1,end-1) = 1;
    AA(end-1,end-3) = -1;
    BB(end-1) = 0;

    TT = AA\BB;
    Tm = TT;

    Tfluido = TT(2:2:2*Nx);
    Tsolido = TT(1:2:2*Nx);
           if abs(time-t1)<=dt/2;
                 iplot=1;
           elseif abs(time-t2)<=dt/2;
                 iplot=1;
           elseif abs(time-t3)<=dt/2;
                 iplot=1; 
            elseif abs(time-t4)<=dt/2;
                 iplot=1; 
           end
           if iplot==1
                 plot(xx,Tsolido,'linewidth',2);
                 hold on;
                 plot(xx,Tfluido,'--','linewidth',2);
                 iplot=0;
                 hold on;
                 
           end

end




