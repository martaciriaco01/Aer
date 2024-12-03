clear all
close all
%% INPUT
body=2; % 0 -> circular cylinder; 1-> rectangular cylinder; 2-> airfoil

plot_fig=1; %1 per plot

U = 1; %Velocità indisturbata
alfa = 2*pi/180; % angolo di incidenza
% Scomposizione della velocità nelle due componenti
Uinf(1)=U*cos(alfa); Uinf(2)=U*sin(alfa);

%% LOAD POINTS: Panel extrema
if body==0 % CIRCULAR CYLINDER:
    Np = 150; R=1; Xc = 0; Yc = 0;
    Xlo = -0.3; Ylo=-0.3;
    Xup = 0.3; Yup = 0.3;
    %Estremi pannelli
    xe = zeros(Np+1,1); ye = xe;
    for i = 0: Np
        xe(i+1) = Xc+R*cos(2*pi/Np*i);
        ye(i+1) = Yc - R*sin(2*pi/Np*i);
    end
elseif body==1     % RECTANGULAR CYLINDER
    Npx=40; Npy=8; Np=2*(Npx+Npy);
    D=1; L=5*D; Xc=0; Yc=0;
    Xup=Xc+L/2; Xlo=Xc-L/2; Yup=Yc+D/2; Ylo=Yc-D/2;
    %Estremi pannelli
    Npy2 = floor(Npy/2);
    xe = zeros(Np+1,1); ye=xe;
    for i = 0 : Npy2
        xe(i+1) = Xup;
        ye(i+1) = Yup-(Yup-Ylo)/2-(Yup-Ylo)/Npy*i; 
    end
    for i = Npy2+1: Npy2+Npx
        xe(i+1) = Xup-(Xup-Xlo)/Npx*(i-Npy2);
        ye(i+1)= Ylo;
    end
    for i = Npy2+Npx+1: Npy2+Npx+Npy
        xe(i+1)=Xlo;
        ye(i+1)=Ylo+(Yup-Ylo)/Npy*(i-Npy2-Npx);
    end
    for i = Npy2+Npx+Npy+1 : Npy2+Npx+Npy+Npx
        xe(i+1) = Xlo+(Xup-Xlo)/Npx*(i-Npy2-Npx-Npy);
        ye(i+1)=Yup;
    end
    for i = Npy2 + Npx + Npy + Npx + 1 : Npy2 + Npx + Npy + Npx + Npy2 
        xe(i+1)=Xup;
        ye(i+1)=Yup-(Yup-Ylo)/Npy*(i-Npy2-2*Npx-Npy);
    end
elseif body==2 %AIRFOIL
    %M = 0; P=0; SS=24; n=100;
    M = 0; P=0; SS=12; n=50;

    % M  : ordinata massima della linea media             (una cifra)
    %
    % P  : posizione sulla corda dell'ordinata massima M  (una cifra)
    %
    % SS : spessore massimo del profilo simmetrico        (due cifre)
    %
    % n  : numero di punti lungo la corda distribuiti uniformemente,
    %      inclusi i punti del bordo di attacco e del bordo di uscita,  
    %      in corrispondenza dei quali si calcolano i punti sul dorso
    %      e sul ventre
    % 
    %
    % M, P e SS possono essere forniti come NUMERI INTERI, compresi, 
    %      rispettivamente,
    %
    %      M  fra 1 e 9,  in  %  della corda      
    %      P  fra 1 e 9   in 1/10  della corda   
    %      SS fra 1 e 99  in  %  della corda
    %
    %      oppure, se si preferisce, come VALORI REALI  < 1  
    %      frazione della corda
    [xd, yd, xv, yv, xe, ye] = NACA_4d_a (M, P, SS, n, 'cos');
    
%    xe=[xe; xd(2:end)];
%    ye=[flipud(yv(:,1)); yd(2:end,1)];
    Np = size(xe,1)-1;

    %figure(1); hold on; box on;
    plot(xe,ye,'r*')
%% 

end
%% 




%% Build geometry information needed
% Punti di Controllo dei pannelli
xp = 0.5*(xe(1:Np) + xe(2:Np+1));
yp = 0.5*(ye(1:Np) + ye(2:Np+1));

% luneghezze dei pannelli
length= sqrt((ye(2:Np+1)-ye(1:Np)).^2+(xe(2:Np+1)-xe(1:Np)).^2);
% normals
normal = zeros(Np,2); % In column 1 we put x component, in column 2 the y component
normal(:,1) = -(ye(2:Np+1)-ye(1:Np))./length;
normal(:,2)=   (xe(2:Np+1)-xe(1:Np))./length; 
% tangent vector
tangent = zeros(Np,2);
tangent(:,1) = (xe(2:Np+1)-xe(1:Np))./length;
tangent(:,2) = (ye(2:Np+1)-ye(1:Np))./length;




%% Build the linear system Ax=b
A = zeros( Np+1, Np+1);
b = zeros( Np+1, 1);

% non penetration boundary condition
for i = 1: Np
    for j = 1: Np
        [sou.u,sou.v] = source(xe(j),ye(j),xe(j+1),ye(j+1),xp(i),yp(i));
        [vor.u,vor.v] = vortex( xe(j),ye(j),xe(j+1),ye(j+1),xp(i),yp(i));  
        A(i,j) = sou.u*normal(i,1) + sou.v*normal(i,2);
        A(i,Np+1) = A(i,Np+1) + vor.u*normal(i,1) + vor.v*normal(i,2);
    end
end

%Kutta condition
for j = 1:Np
    [sou.u,sou.v] = source(xe(j),ye(j),xe(j+1),ye(j+1),xp(1),yp(1)); 
    [vor.u,vor.v] = vortex(xe(j),ye(j),xe(j+1),ye(j+1),xp(1),yp(1));
    A(Np+1,j) = sou.u*tangent(1,1)+sou.v*tangent(1,2);
    A(Np+1,Np+1) = A(Np+1,Np+1) + vor.u*tangent(1,1) + vor.v*tangent(1,2);

    [sou.u,sou.v] = source(xe(j),ye(j),xe(j+1),ye(j+1),xp(Np),yp(Np));    
    [vor.u,vor.v] = vortex(xe(j),ye(j),xe(j+1),ye(j+1),xp(Np),yp(Np)); 
    A(Np+1,j) = A(Np+1,j) +  sou.u*tangent(Np,1)+sou.v*tangent(Np,2);
    A(Np+1,Np+1) = A(Np+1,Np+1) + vor.u*tangent(Np,1) + vor.v*tangent(Np,2);    
end

% right hand side
for i = 1:Np
    b(i) = -Uinf(1)*normal(i,1) - Uinf(2)*normal(i,2);
end
b(Np+1) = -Uinf(1)*tangent(1,  1) - Uinf(2)*tangent(1,  2) ...
                 -Uinf(1)*tangent(Np,1) - Uinf(2)*tangent(Np,2);

%% Solve the linear system
sol = linsolve(A,b);

sigma = sol(1:Np);
gamma= sol(Np+1);


%% Compute velocity on control points and pressure distribution
velu= zeros(Np,1); velv = velu;
for i =1:Np
    velu(i) = Uinf(1); velv(i)=Uinf(2);
    for j = 1:Np
        [sou.u,sou.v] = source(xe(j),ye(j),xe(j+1),ye(j+1),xp(i),yp(i));
        [vor.u,vor.v] = vortex(xe(j),ye(j),xe(j+1),ye(j+1),xp(i),yp(i));
       velu(i) = velu(i) + sigma(j)*sou.u + gamma*vor.u;
       velv(i) = velv(i) + sigma(j)*sou.v + gamma*vor.v;
    end
end

if (max(velu.*normal(:,1) + velv.*normal(:,2))>10^(-14))
    disp('There is a bug in the program!')
    stop
end

Vt = velu.*tangent(:,1) + velv.*tangent(:,2);
Vn = velu.*normal(:,1) + velv.*normal(:,2);

Cp = 1-Vt.^2/U^2;

% Cl = sum (Cp_i * L_i* cos(theta_i) )
Cl = -Cp'*( length.*normal(:,2) )

circ= sum(length.*gamma)
rho=1;
Lift = rho * U* circ
Cl1 = Lift/(0.5*rho*U^2)






%% plot velocity field

if plot_fig==1
    
    figure(1); hold on; box on;
    plot(xe,ye,'r*')
    plot(xe,ye,'k')
    plot(xp,yp,'bs')
    axis equal
    
    
    Nfx=100; Nfy = 100;
    x=linspace(-4,4,Nfx);
    y=linspace(-1.5,1.5,Nfy);

    Ufield= zeros(Nfx,Nfy); Vfield = Ufield;

    if (body==2)
    ydplus = interp1(xd(:,1),yd(:,1),x);
    yvplus = interp1(xv(:,1),yv(:,1),x);
    end
    for i =1:Nfx
        for j = 1:Nfy
            if (body==0 && ((x(i)-Xc)^2+(y(j)-Yc)^2>R^2)) 
                Ufield(i,j) = Uinf(1); Vfield(i,j) = Uinf(2);
               for ip = 1:Np
                  [sou.u,sou.v] = source(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  [vor.u,vor.v] = vortex(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  Ufield(i,j) = Ufield(i,j) + sigma(ip)*sou.u + gamma*vor.u;
                  Vfield(i,j) = Vfield(i,j) + sigma(ip)*sou.v + gamma*vor.v;
               end
            elseif (body==1 && (x(i)< Xlo || y(j)<Ylo || x(i)>Xup|| y(j)>Yup))
                Ufield(i,j) = Uinf(1); Vfield(i,j) = Uinf(2);
               for ip = 1:Np
                  [sou.u,sou.v] = source(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  [vor.u,vor.v] = vortex(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  Ufield(i,j) = Ufield(i,j) + sigma(ip)*sou.u + gamma*vor.u;
                  Vfield(i,j) = Vfield(i,j) + sigma(ip)*sou.v + gamma*vor.v;
               end
            elseif (body==2 && ( x(i)<xd(1,1) || x(i)>xd(end,1) || ( ( x(i)<xd(end,1) && x(i)>xd(1,1) ) && ( y(j)>ydplus(i) || y(j)<yvplus(i)  )) ) )
                Ufield(i,j) = Uinf(1); Vfield(i,j) = Uinf(2);
               for ip = 1:Np
                  [sou.u,sou.v] = source(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  [vor.u,vor.v] = vortex(xe(ip),ye(ip),xe(ip+1),ye(ip+1),x(i),y(j));
                  Ufield(i,j) = Ufield(i,j) + sigma(ip)*sou.u + gamma*vor.u;
                  Vfield(i,j) = Vfield(i,j) + sigma(ip)*sou.v + gamma*vor.v;
               end
            end            
        end
    end
    [X,Y] = meshgrid(x,y);
    figure(3); hold on; box on;
    contourf(X',Y',Ufield,100,'LineStyle','None');
    colormap(flipud(hot));
    plot(xe,ye,'Linewidth',2,'Color','k')
    colorbar('off')
    view(2);
    tx=xlabel('$x$');
    ty=ylabel('$y$');
    tx.Interpreter='latex';
    ty.Interpreter='latex';
    set(gca,'TickLabelInterpreter', 'latex');  
    x0=10;
    y0=10;
    width=1250; 
    height=580;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca;
    axis equal
    axis off
    ax.FontSize = 40;

    figure(4); hold on; box on;
    contourf(X',Y',Vfield,100,'LineStyle','None');
    colormap(bluewhitered(256)); colorbar;
    c=colorbar('northoutside');
    c.TickLabelInterpreter='latex';
    colorbar('off')
    plot(xe,ye,'Linewidth',2,'Color','k')
    view(2);
    tx=xlabel('$x$');
    ty=ylabel('$y$');
    tx.Interpreter='latex';
    ty.Interpreter='latex';
    set(gca,'TickLabelInterpreter', 'latex');  
    x0=10;
    y0=10;
    width=1250; 
    height=580;
    set(gcf,'position',[x0,y0,width,height])
    ax = gca;
    axis equal
    axis off
    ax.FontSize = 40;
    
    if body==2
        figure(); hold on; box on; grid on
        plot(xp,-Cp,'LineWidth',1.5,'Color','r')
        tx=xlabel('$x$');
        ty=ylabel('$-Cp$');
        tx.Interpreter='latex';
        ty.Interpreter='latex';
        set(gca,'TickLabelInterpreter', 'latex');  
        x0=10;
        y0=10;
       width=1250; 
       height=580;
       set(gcf,'position',[x0,y0,width,height])
       ax = gca;
       ax.FontSize = 40;
    end
end

%% Source and vortex functions
function [u,v] = source( xe1, ye1, xe2, ye2, x, y)
l=sqrt((xe2-xe1)^2+(ye2-ye1)^2); %lunghezza del pannello
% mi metto nel riferimento locale
xl2 = l; yl2=0;
xl=(x-xe1)*(xe2-xe1)/l+(y-ye1)*(ye2-ye1)/l; yl=-(x-xe1)*(ye2-ye1)/l+(y-ye1)*(xe2-xe1)/l; %posizione del punto (x,y) nel riferimento locale
r1=sqrt(xl^2+yl^2); r2=sqrt((xl-xl2)^2+(yl-yl2)^2); %calcolo r1 e r2
theta1=atan2(yl,xl); theta2=atan2(yl-yl2,xl-xl2); % calcolo theta1 e theta2
% attenzione ad errori numerici che possono far sbagliare la determinazione
% dell'angolo.
if (abs(theta1)<10^(-12) && abs(theta2)>3); theta1=0; theta2=pi; end
if (abs(theta2)<10^(-12) && abs(theta1)>3); theta2=0; theta1=-pi; end
% velocità indotta nel riferimento locale
vl=1/(2*pi)*(theta2-theta1); ul=1/(2*pi)*log(r1/r2);
% velocità indotta nel riferimento globale
u=ul*(xe2-xe1)/l-vl*(ye2-ye1)/l; v=ul*(ye2-ye1)/l+vl*(xe2-xe1)/l;
end

function [u,v] = vortex( xe1, ye1, xe2, ye2, x, y)
l=sqrt((xe2-xe1)^2+(ye2-ye1)^2);
xl2 = l; yl2=0;
xl=(x-xe1)*(xe2-xe1)/l+(y-ye1)*(ye2-ye1)/l; yl=-(x-xe1)*(ye2-ye1)/l+(y-ye1)*(xe2-xe1)/l;
r1=sqrt(xl^2+yl^2); r2=sqrt((xl-xl2)^2+(yl-yl2)^2);
theta1=atan2(yl,xl); theta2=atan2(yl-yl2,xl-xl2);
if (abs(theta1)<10^(-12) && abs(theta2)>3); theta1=0; theta2=pi; end
if (abs(theta2)<10^(-12) && abs(theta1)>3); theta2=0; theta1=-pi; end
ul=1/(2*pi)*(theta2-theta1); vl=-1/(2*pi)*log(r1/r2);
u=ul*(xe2-xe1)/l-vl*(ye2-ye1)/l; v=ul*(ye2-ye1)/l+vl*(xe2-xe1)/l;
end


  


