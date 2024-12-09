%% traccia Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 1;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha));
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y];
U_inf_normal = [-U_inf(2); U_inf(1)];
U_inf_normal = U_inf_normal ./ norm(U_inf_normal);

TestCase = 0;

CodiceProfilo = '0012';
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
% [x,y]=createProfile(CodiceProfilo,NPannelli,Chord);

Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat'));
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

figure;
plot(x, y, 'o-')
axis equal

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As

for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)';
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna

            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));


%% Risoluzione sistema lineare

Soluzione = linsolve(matriceA,TermineNoto); % cj (101) e cj+1 (1) --> 102


%% Calcolo del cp e della velocità sui pannelli

velu= zeros(NPannelli,1); 
velv = velu;

for i = 1 : NPannelli
    velu(i) = U_inf_x; 
    velv(i) = U_inf_y;
    Centropan = Centro(i,:)';
    for j = 1 : NPannelli
        Estremo1pan = Estremo_1(j,:)';
        Estremo2pan = Estremo_2(j,:)';
        L2G_TransfMatrix_pan = squeeze(L2G_TransfMatrix(j, :, :));
        G2L_TransfMatrix_pan = squeeze(G2L_TransfMatrix(j, :, :));
        [U_s] = ViSorgente(Centropan, Estremo1pan, Estremo2pan, L2G_TransfMatrix_pan, G2L_TransfMatrix_pan);
        [U_v] = ViVortice(Centropan, Estremo1pan, Estremo2pan, L2G_TransfMatrix_pan, G2L_TransfMatrix_pan);
       velu(i) = velu(i) + Soluzione(j)*U_s(1) + Soluzione(NPannelli+1)*U_v(1);
       velv(i) = velv(i) + Soluzione(j)*U_s(2) + Soluzione(NPannelli+1)*U_v(2);
    end
end

if (max(velu.*Normale(:,1) + velv.*Normale(:,2))>10^(-14))
    disp('There is a bug in the program!')
    return
end

Vt = velu.*Tangente(:,1) + velv.*Tangente(:,2);
Vn = velu.*Normale(:,1) + velv.*Normale(:,2);

U_inf = 1; % velocità indisturbata

Cp = 1-Vt.^2/U_inf^2;



%% COEFFICIENTE DI PORTANZA --> Cl = sum (Cp_i * L_i* cos(theta_i) )

% coefficiente di portanza mediante il coefficiente di pressione
Cl = -Cp' * ( lunghezza' .* Normale(:,2) );

% coefficiente di portanza mediante circolazione
circ= sum(lunghezza' .* Soluzione(NPannelli+1)); % circolazione totale
rho = 1;
Lift = rho * U_inf * circ;  % teorema di Jutta-Joukowsky per calcolo portanza
Cl1 = Lift/(0.5*rho*U_inf^2);



%% Cm rispetto al LE

% Coordinate dei centri dei pannelli
x_centro = Centro(:, 1); 
z_centro = Centro(:, 2); 

r_cross_n = x_centro .* Normale(:,2) - z_centro .* Normale(:,1); 
Cm_LE = Cp' * (lunghezza' .* r_cross_n);


%% Cm rispetto al CENTRO AERODINAMICO

% Coordinate dei centri dei pannelli
x_centro = Centro(:, 1); 
z_centro = Centro(:, 2); 

% Centro aerodinamico a x = 1/4 
x_c4 = x_centro - 1/4; 
 
r_c4_cross_n = x_c4 .* Normale(:,2) - z_centro .* Normale(:,1);
Cm_c4 = (Cp' * (lunghezza' .* r_c4_cross_n));



%% creo un vettore con Cl,Cl1,Cm_c4 per vedere subito i risultati

risultati = [Cl, Cl1, Cm_c4]





