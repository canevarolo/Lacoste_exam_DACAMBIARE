% Esercitazione 14, esercizio 3
% Simone Canevarolo
% S269893
% 20/02/2024

clear all
close all
clc

ss1 = 2e-3; % spessore lastra 1, m
ss2 = 3e-3; % spessore lastra 2, m
ssliq = 5e-3; % spessore liquido, m
ll = 2; % lunghezza condotto, m

As = 
Vol = 

roacc = 7800; % densità, kg/m^3
cpacc = 500; % calore specifico acciaio, J/Kg/K
kk = 30; % condicibilità, W/m/K
hh = 1e3; % coefficiente di scambio termico, W/m^2/K

vv = 1e-2; % velocità, m/s

Teq = 25+273; % temperatura di equilibrio, K
Tinh2o = Teq; % temperatura di ingresso dell'acqua nel condotto

q0 = 3e3; % 

dx = 1e-3; % m
xx = (0:dx:ll)';
Nx = length(xx);

qvol = @(x) q0*exp(-5*(x-1).^2);

figure(1)
plot(xx,qvol(xx),'linewidth',2)

%%

a1 = kk*dt/roacc/cpacc/dx^2;
b1 = hh*dt*As/roacc/cpacc/Vol;
b2 =  hh*dt*As/roacc/cpacc/Vol;
bw1 =  hh*dt*As/roacc/cpacc/Vol;
bw2 =  hh*dt*As/roacc/cpacc/Vol;
cc = vv*dt/dx;
q1 = qvol(xx)*dt/roacc/cpacc;

% Assemblaggio matrice

supsupsupdiag = zeros(3*Nx,1);
supsupdiag = zeros(3*Nx,1);
supdiag = zeros(3*Nx,1);
maindiag = zeros(3*Nx,1);
subdiag = zeros(3*Nx,1);
subsubdiag = zeros(3*Nx,1);
subsubsubdiag = zeros(3*Nx,1);

supsupsupdiag(1:3:end) = -a1;
supsupsupdiag(2:3:end) = -cc;
supsupsupdiag(3:3:end) = -a1;

supdiag(2:3:end) = -bw1;
supdiag(3:3:end) = -b2;

maindiag(1:3:end) = 1+2*a1+b1;
maindiag(2:3:end) = 1+cc+bw1+bw2;
maindiag(3:3:end) = 1+2*a1+b2;

subdiag(1:3:end) = -b1;
subdiag(2:3:end) = -bw2;

subsubsubdiag(1:3:end) = -a1;
subsubsubdiag(3:3:end) = -a1;

Band = [subsubsubdiag,subsubdiag,subdiag,maindiag,supdiag,supsupdiag,supsupsupdiag];

AA = spdiags(Band,-3:3,3*Nx,3*Nx);

dt = 1;
Tguess = Tinh2o;

bb1 = -q1;
bb2 = zeros(Nx,1);
bb3 = zeros(Nx,1);

bb = zeros(3*Nx,1);
bb(1:3:end) = bb1;
bb(2:3:end) = bb2;
bb(3:3:end) = bb3;


for ii = 1:1000

    % ii = ii+1;
    time = time+dt;

    AA

    

end




