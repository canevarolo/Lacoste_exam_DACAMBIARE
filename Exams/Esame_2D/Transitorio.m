clear all
clear vars
clc

%Dati
Ly=25e-2;
Lx=15e-2;
Lx1=Lx/3;
Lx2=Lx-Lx1;
Ta=20;
ha=20;
Tc=5;
cond=2;

%% Mappa di temperatura in condizioni stazionarie
%Discretizzazione
dx=5e-3;
dy=dx;

x1=(0:dx:Lx1)';
nx1=length(x1);
x2=(0:dx:Lx2)';
nx2=length(x2);

xx=(0:dx:Lx)';
nx=length(xx);
yy=(0:dy:Ly)';
ny=length(yy);

ntot=nx*ny;

[xmat,ymat]=meshgrid(xx,yy);

%Matrice
AA=sparse([],[],[],ntot,ntot,5*ntot);
bb=zeros(ntot,1);

for ii=2:nx-1
    for jj=2:ny-1
        kk=ii+(jj-1)*nx;
        kn=kk-1;
        ks=kk+1;
        ko=kk-nx;
        ke=kk+nx;
        AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
        AA(kk,kn)=cond/dx^2;
        AA(kk,ks)=cond/dx^2;
        AA(kk,ko)=cond/dy^2;
        AA(kk,ke)=cond/dy^2;
    end
end

%Bordi
%Nord condizione di neumann adiabatica
% kk=1:nx:ntot;
% AA(kk,:)=0;
ii=1;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,ks)=2*cond/dx^2;
    AA(kk,ke)=cond/dy^2;
    AA(kk,ko)=cond/dy^2;
    bb(kk)=0;
end

%Sud condizione di Robin
ii=nx;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2)-ha/dx;
    AA(kk,kn)=2*cond/dx^2;
    AA(kk,ko)=cond/dy^2;
    AA(kk,ke)=cond/dy^2;
    bb(kk)=-ha/dx*Ta;
end

%Ovest Neumann adiabatico
jj=1;
for ii=2:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,ks)=cond/dx^2;
    AA(kk,kn)=cond/dx^2;
    AA(kk,ke)=2*cond/dy^2;
    bb(kk)=0;
end

%Est divisa tra adiabatica e temperatura imposta
jj=ny;
for ii=2:nx1-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,kn)=cond/dx^2;
    AA(kk,ks)=cond/dx^2;
    AA(kk,ko)=2*cond/dy^2;
    bb(kk)=0;
end
for ii=nx1:nx2
    kk=ii+(jj-1)*nx;
    AA(kk,kk)=1;
    bb(kk)=Tc;
end
for ii=nx2+1:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    AA(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    AA(kk,kn)=cond/dx^2;
    AA(kk,ks)=cond/dx^2;
    AA(kk,ko)=2*cond/dy^2;
    bb(kk)=0;
end

%Nodi
kk=nx;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
AA(kk,kk-1)=cond/dx^2;
AA(kk,kk+nx)=cond/dy^2;
bb(kk)=-ha/dx*Ta;
kk=ntot;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
AA(kk,kk-1)=cond/dx^2;
AA(kk,kk-nx)=cond/dy^2;
bb(kk)=-ha/dx*Ta;
kk=1;
AA(kk,:)=0;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2);
AA(kk,kk+1)=cond/dx^2;
AA(kk,kk+nx)=cond/dy^2;
bb(kk)=0;
kk=ntot-nx+1;
AA(kk,kk)=-cond*(1/dx^2+1/dy^2);
AA(kk,kk+1)=cond/dx^2;
AA(kk,kk-nx)=cond/dy^2;
bb(kk)=0;

TT=AA\bb;

%% Transitorio
dt=1800;
rho=8000;
cp=500;
BB=eye(length(AA))-AA*dt/rho/cp;

%Bordi
%Nord condizione di neumann adiabatica
% kk=1:nx:ntot;
% AA(kk,:)=0;
ii=1;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    BB(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    BB(kk,ks)=2*cond/dx^2;
    BB(kk,ke)=cond/dy^2;
    BB(kk,ko)=cond/dy^2;
end

%Sud condizione di Robin
ii=nx;
for jj=2:ny-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    BB(kk,kk)=-2*cond*(1/dx^2+1/dy^2)-ha/dx;
    BB(kk,kn)=2*cond/dx^2;
    BB(kk,ko)=cond/dy^2;
    BB(kk,ke)=cond/dy^2;
end

%Ovest Neumann adiabatico
jj=1;
for ii=2:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    BB(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    BB(kk,ks)=cond/dx^2;
    BB(kk,kn)=cond/dx^2;
    BB(kk,ke)=2*cond/dy^2;
end

%Est divisa tra adiabatica e temperatura imposta
jj=ny;
for ii=2:nx-1
    kk=ii+(jj-1)*nx;
    kn=kk-1;
    ks=kk+1;
    ko=kk-nx;
    ke=kk+nx;
    BB(kk,kk)=-2*cond*(1/dx^2+1/dy^2);
    BB(kk,kn)=cond/dx^2;
    BB(kk,ks)=cond/dx^2;
    BB(kk,ko)=2*cond/dy^2;
end
%Nodi
kk=nx;
BB(kk,:)=0;
BB(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
BB(kk,kk-1)=cond/dx^2;
BB(kk,kk+nx)=cond/dy^2;

kk=ntot;
BB(kk,:)=0;
BB(kk,kk)=-cond*(1/dx^2+1/dy^2)-ha/dx;
BB(kk,kk-1)=cond/dx^2;
BB(kk,kk-nx)=cond/dy^2;

kk=1;
BB(kk,:)=0;
BB(kk,kk)=-cond*(1/dx^2+1/dy^2);
BB(kk,kk+1)=cond/dx^2;
BB(kk,kk+nx)=cond/dy^2;

kk=ntot-nx+1;
BB(kk,:)=0;
BB(kk,kk)=-cond*(1/dx^2+1/dy^2);
BB(kk,kk+1)=cond/dx^2;
BB(kk,kk-nx)=cond/dy^2;

Told=TT;
tt=0:dt:3600*24*3;
for gg=1:length(tt)
    dd=Told+bb*dt/(rho*cp);
    %Nord condizione di neumann adiabatica
    % kk=1:nx:ntot;
    % AA(kk,:)=0;
    ii=1;
    for jj=2:ny-1
        kk=ii+(jj-1)*nx;
        dd(kk)=0;
    end
    
    %Sud condizione di Robin
    ii=nx;
    for jj=2:ny-1
        kk=ii+(jj-1)*nx;
        dd(kk)=-ha/dx*Ta;
    end
    
    %Ovest Neumann adiabatico
    jj=1;
    for ii=2:nx-1
        kk=ii+(jj-1)*nx;
        dd(kk)=0;
    end
    
    %Est divisa tra adiabatica e temperatura imposta
    jj=ny;
    for ii=2:nx-1
        kk=ii+(jj-1)*nx;
        dd(kk)=0;
    end
    %Nodi
    kk=nx;
    dd(kk)=-ha/dx*Ta;
    kk=ntot;
    dd(kk)=-ha/dx*Ta;
    kk=1;
    dd(kk)=0;
    kk=ntot-nx+1;
    dd(kk)=0;
    T0=BB\dd;
    Told=T0;
end
figure(1)
Tplot=reshape(T0,nx,ny);
surf(xmat',ymat',Tplot)
xlabel('X(m)')
ylabel('Y(m)')
zlabel('T (°C)')