% Esercitazione 4, esercizio 6
% Simone Canevarolo
% S269893
% 11/11/2023

clear all
close all
clc

%%

kma = 0.55; % W/(m*K) conduttività termica mattone
kis = 0.04; % W/(m*K) conduttività termica isolante

llma = 25e-2; % m spessore mattoni
llis = 10e-2; % m spessore isolante

Test = 268; % K temperatura esterna
hest = 25; % W/(m*K) coefficiente di scambio termico esterno
Tin = 293; % K temperatura interna

%% Parto con il caso di "omogeneità"

dx = 1e-4; % m
xx = (0:dx:(llma+llis))';
NN = length(xx);

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = ones(NN,1);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];

AA = spdiags(Band,-1:1,NN,NN);

bb = zeros(NN,1);

% Condizioni al contorno

% Robin
AA(1,1) = -1-hest*dx/kis;
AA(1,2) = 1;
bb(1) = -Test*hest*dx/kis;

% Dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tin;

TT = AA\bb;

%% Affronto ora il caso di interfaccia

kma = 0.55; % W/(m*K) conduttività termica mattone
kis = 0.04; % W/(m*K) conduttività termica isolante

llma = 25e-2; % m spessore mattoni
llis = 10e-2; % m spessore isolante

Test = 268; % K temperatura esterna
hest = 25; % W/(m*K) coefficiente di scambio termico esterno
Tin = 293; % K temperatura interna

xxis = (0:dx:llis)';
xxma = (llis+dx:dx:llma)';

NNis = length(xxis);
% NNma = length(xxma);

sub_diag = ones(NN,1);
main_diag = -2*ones(NN,1);
sup_diag = ones(NN,1);

Band = [[sub_diag(2:end);0],main_diag,[0;sup_diag(1:end-1)]];
AA = spdiags(Band,-1:1,NN,NN);
bb = zeros(NN,1);

% Robin a x=0
AA(1,1) = -1-hest*dx/kis;
AA(1,2) = 1;
bb(1) = -Test*hest*dx/kis;

% Robin all'interfaccia
AA(NNis,NNis-1) = -kis;
AA(NNis,NNis) = kis+kma;
AA(NNis,NNis+1) = -kma;
bb(NNis) = 0;

% Dirichlet
AA(end,end-1) = 0;
AA(end,end) = 1;
bb(end) = Tin;

TT = AA\bb;

Tmax = max(TT)
Tint = TT(NNis)

%% Grafico

figure(1)
plot(xx,TT)

title('Andamento della temperatura con discontinuità')
xlabel('Lunghezza [m]')
ylabel('Temperatura [K]')