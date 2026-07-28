% Tema d'esame giugno 2023
% Simone Canevarolo
% S269893
% 1/02/2024

clear all
close all
clc

% Trascuro dimensioni spaziali, considero il problema monodimensionale
ll = 0.3; % m

Tsx = 77; % Temperatura sinistra, K
Tsc = 20; % temperatura superconduttore, K
Tmedia = (Tsx+Tsc)/2; % temperatura a ll/2 considerando l'andamento di temperatura lineare nel pezzo, K
epsilon = 0.7; % emissività

dx = 1e-4;
xx = (0:dx:ll)';
Nx = length(xx);

% Metallo 1

k1 = @(T) 1/3*(-1.15e-6*T^4 +1.18e-3*T^3 -2.62e-1*T^2 +1.83e1*T -38.9);

Tguess = Tmedia*ones(Nx,1);

toll = 1e-5;
err = toll+1;
ii = 1;
maxii = 1000;

while err>toll && ii<maxii

        ii = ii+1;

        



end


