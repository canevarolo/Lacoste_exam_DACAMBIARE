% Esame 28 gennaio 2022
% Simone Canevarolo
% S269893
% 2/07/2024

clear all
close all
clc

lungh = 3; % m
raggio = 1e-2;
ss = 3e-3;

port = 0.02; % Kg/s
rovol = 1.4e-3; % Kg/m^3
cppvc = 1400; % J/Kg/K
kpvc = 0.16; % W/m/K

Taria = 278; % K
hout = 10; % W/m^2/K

Tin = 50; % K

roh2o = 1e3; % kg/m^3
Abase = pi*(raggio^2);

uu = port*Abase/roh2o; % m/s
Asint = 2*pi*raggio*lungh; % m^2
Asest = 2*pi*(raggio+ss)*lungh; % m^2
Volint = Abase*lungh; % m^3
Volest = pi*(raggio+ss)^2*lungh-Volint; % m^3

dx = 1e-3; % m
xx1 = (0:dx:lungh)';
xx2 = xx1;
Nx1 = length(xx1);

xx = [xx1;xx2];
Nx = length(xx);

aa = rovol*cppvc*uu/dx;
cc = hout*Asint/Volint;

supsup = zeros(Nx,1);
sup = zeros(Nx,1);
main = zeros(Nx,1);
sub = zeros(Nx,1);
subsub = zeros(Nx,1);

supsup(1:2:end) = kpvc/dx^2*ones(Nx1,1);
supsup(2:2:end) = zeros(Nx1,1);

sup(1:2:end,1) = zeros(Nx1,1);
sup(2:2:end) = zeros(Nx1,1);

main(1:2:end,1) = ((-2*kpvc/dx^2+aa)+hout*Asest/Volest)*ones(Nx1,1);
main(2:2:end,1) = (-aa-cc)*ones(Nx1,1);

sub(1:2:end,1) = zeros(Nx1,1);
sub(2:2:end) = -hout*Asint/Volint;

subsub(1:2:end) = kpvc/dx^2*ones(Nx1,1);
subsub(2:2:end,1) = aa*ones(Nx1,1);

Band = [subsub, sub, main, sup, supsup];

AA = spdiags(Band,-2:2,Nx,Nx);
bb = zeros(Nx,1);
bb(1:2:end) = Taria*hout*Asest/Volest;

% Neumann adiabatico sx
AA(1,:) = 0;
AA(1,1) = 1;
AA(1,3) = -1;
bb(1) = 0;

% Dirichlet dx solido

% AA() = 1;
% bb(end) = Tin;

% Dirichlet dx liquido
AA(end,:) = 0;
AA(end,end) = 1;
bb(end) = Tin;

TT = AA\bb;

TT1 = TT(1:2:end);
TT2 = TT(2:2:end);

figure(1)
plot(xx1,TT1)
hold on
plot(xx1,TT2)



