% Esercitazione 7, esercizio 1
% Simone Canevarolo
% S269893
% 27/11/2023

clear all
close all
clc

%%

ll = 50e-3; % lunghezza m
kk = 5; % conducibilità termica W/(m*K)
Ta = 300; % temperatura imposta K
Tb = 300; % temperatura imposta K

qv = @(x) 5e5*sin(pi*x/ll);

dx = 1e-5;
xx = (0:dx:ll)';
NN = length(xx);

%%

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = sub_diag;

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = -qv(xx)*dx^2/kk;

%%

% Dirichlet Ta

AA(1,1) = 1;
AA(1,2) = 0;
bb(1) = Ta;

% Dirichlet Tb

AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tb;

TT = AA\bb;

%%

figure(1)
plot(xx,TT,'linewidth',2)
xlabel('lunghezza [m]')
ylabel('temperatura [K]')