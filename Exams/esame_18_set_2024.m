% Considerando sempre perfetto l’accoppiamento termico con l’azoto liquido, e il coefficiente di scambio 
% termico con il vapore di azoto pari a 50 × (𝑇𝐻𝑂𝑇𝑆𝑃𝑂𝑇(𝑡) − 70) W/m2K, scrivere un report dettagliato 
% sull’evoluzione del transitorio termico della sbarra fino al raggiungimento di un nuovo stato stazionario. 
% In particolare, dopo aver correttamente formulato il problema matematico di riferimento e discusso gli 
% schemi numerici utilizzati per la discretizzazione in spazio e tempo dell’equazione di riferimento e relative 
% condizioni iniziali e al contorno, si studi per la sezione 2D mostrata in figura:

% a. La temperatura massima 𝑇𝐻𝑂𝑇𝑆𝑃𝑂𝑇 all’interno della sbarra, producendo e commentando un grafico 
% che mostri l’evoluzione di 𝑇𝐻𝑂𝑇𝑆𝑃𝑂𝑇 nel tempo.

% b. La distribuzione di temperatura nella sbarra, producendo la mappa di temperatura nello stazionario 
% finale e commentando se e come rispetti le condizioni al contorno del problema.

% c. La distribuzione di temperatura lungo l’asse di simmetria verticale della sbarra, producendo e 
% commentando un grafico che mostri contemporaneamente il profilo di temperatura ivi calcolato 
% dopo 1s, 2 s, 5 s e 10 s dall’inizio del transitorio.

% d. L’accuratezza con cui viene calcolato il valore di 𝑇𝐻𝑂𝑇𝑆𝑃𝑂𝑇 nel nuovo stazionario, al variare dei 
% parametri di discretizzazione spaziale, discutendo se e come possa essere raggiunta una 
% accuratezza di 0.1 K.

% e. Il bilancio di potenza nello stazionario finale, verificando e commentando come la potenza 
% istantanea che la sbarra trasferisce all’azoto liquido e al suo vapore eguagli la potenza generata per 
% effetto Joule all’interno della sbarra.


% Esame 18/09/2024
% Simone Canevarolo
% s269893

clear all
close all
clc

rovol = 8e3;
cp = 120; 
cond = 15;
roel = 8e-7;
Tvap = 77;
II = 350; % A

base = 1e-2;
altezza = 2e-2;
llvap = 1e-2;

dx = 1e-3;
dy = dx;

xx = (0:dx:altezza)';
yy = (0:dy:base)';
xvap = (0:dx:llvap)';

Nx = length(xx);
Ny = length(yy);
Nvap = length(xvap);
Ntot = Nx*Ny;

hhfunz = @(Thot) 50*(Thot-70);
Asez = base*altezza;

qvol = II^2*roel/Asez^2;

[xmat,ymat] = meshgrid(xx,yy);
AA = sparse([],[],[],Ntot,Ntot,5*Ntot);

tempomax = 10;
dt = 0.1;
tt = (0:dt:tempomax);
Nt = length(tt);

Thotspot = Tvap; % considero il perfetto equilibrio termico
hh = hhfunz(Thotspot);

Thotspotvett = Tvap*ones(Nt,1);

tplot = [1 2 5 10];
iplot = round(tplot/dt)+1;
pp = 1;

aa = cond*dt/rovol/cp;

    for ii = 2:Nx-1
        for jj = 2:Ny-1
    
            kk = Nx*(jj-1)+ii;
            AA(kk,kk-Nx) = -aa/dy^2;
            AA(kk,kk-1) = -aa/dx^2;
            AA(kk,kk) = 1+2*aa*(1/dx^2+1/dy^2);
            AA(kk,kk+1) = -aa/dx^2;
            AA(kk,kk+Nx) = -aa/dy^2;

    
        end
    end
    
    Tm = Tvap*ones(Ntot,1);
    bb = qvol*dt/rovol/cp*ones(Ntot,1);


for zz = 2:Nt

    hh = hhfunz(Thotspot);

    BB = Tm+bb;
    
    % Sud Ovest - Dirichlet
    % jj = 1 
    kk = Nvap:Nx;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    BB(kk) = Tvap;
    
    
    % Sud - Dirichlet
    % ii = Nx;
    kk = Nx:Nx:Ntot;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    BB(kk) = Tvap;
    
    % Sud Est - Dirichlet
    % jj = Ny;
    kk = (Ny-1)*Nx+Nvap:Ntot;
    AA(kk,:) = 0;
    AA(kk,kk) = eye(length(kk));
    BB(kk) = Tvap;

    % Nord Ovest - robin
    jj = 1;
    for ii = 2:Nvap-1
    
        kk = ii;
        AA(kk,kk-1) = 1/2/dx^2;
        AA(kk,kk) = -1/dx^2-1/dy^2-hh/dy/cond;
        AA(kk,kk+1) = 1/2/dx^2;
        AA(kk,kk+Nx) = 1/dy^2;
        BB(kk) = -Tvap*hh/dy/cond;
    
    end
    
    % Nord Est - Robin
    jj = Ny;
    for ii = 2:Nvap-1

        kk = Nx*(jj-1)+ii;
        AA(kk,kk-Nx) = 1/dy^2;
        AA(kk,kk-1) = 1/2/dx^2;
        AA(kk,kk) = -1/dx^2-1/dy^2-hh/dy/cond;
        AA(kk,kk+1) = 1/2/dx^2;
        BB(kk) = -Tvap*hh/dy/cond;

    end
    
    % Nord - Robin
    ii = 1;
    for jj = 2:Ny-1

        kk = ii+Nx*(jj-1);
        AA(kk,kk-Nx) = 1/2/dy^2;
        AA(kk,kk) = -1/dx^2-1/dy^2-hh/dx/cond;
        AA(kk,kk+1) = 1/dx^2;
        AA(kk,kk+Nx) = 1/2/dy^2;
        BB(kk) = -Tvap*hh/dx/cond;

    end
    
    % Vertice Nord Ovest
    kk = 1;
    AA(kk,kk) = -1/2/dx^2-1/2/dy^2-hh/dx/cond/2-hh/dy/cond/2;
    AA(kk,kk+1) = 1/2/dx^2;
    AA(kk,kk+Nx) = 1/2/dy^2;
    BB(kk) = -Tvap*hh/cond/dx/2-Tvap*hh/cond/dy/2;

    
    % Vertice Nord Est
    kk = (Ny-1)*Nx+1;
    AA(kk,kk-Nx) = 1/2/dy^2;
    AA(kk,kk) = -1/2/dx^2-1/2/dy^2-hh/dx/2/cond-hh/dy/2/cond;
    AA(kk,kk+1) = 1/2/dx^2;
    BB(kk) = -Tvap*hh/cond/dx/2-Tvap*hh/cond/dy/2;
    
    
    TT = AA\BB;

    TTT = reshape(TT,Nx,Ny);

    Thotspot = max(TT);
    Thotspotvett(zz) = Thotspot;
    Tm = TT;

    Tcentro = TTT(:,round(Ny/2));
    % Tcentro = TTT(:,6);

    if zz == iplot(pp)
    
        figure(3)
        plot(xx,Tcentro,'linewidth',2)
        hold on
        title('Plot temperature linea centrale a vari t')
        xlabel('lunghezza [m]')
        ylabel('Temperatura [K]')
        legend('1 s','2 s','5 s','10 s')

        pp = pp+1;
    
    end

end

figure(1)
plot(tt,Thotspotvett,'linewidth',2)
xlabel('tempo [s]')
ylabel('Temperatura [K]')
title('Andamento temperatura massima sulla superficie')

figure(2)
surf(xmat',ymat',TTT);
title('Andamento stazionario finale')

