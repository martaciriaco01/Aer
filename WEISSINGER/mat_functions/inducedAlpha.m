function [alphaIndw, alphaIndt] = inducedAlpha(U, wing, GAMMAw, tail, GAMMAt)
%  INDUCEDALPHA - Ballielo, Ciriaco
%  Restituisce il Vettore Contenente il Valore dell'Angolo
%  di Incidenza Indotto Valutato a 1/4 della Corda dei Pannelli
%  dell'Ultimo Settore Trasversale della Superficie Alare,
%  Dovuto ai Vortici Semi-Infiniti della Scia.
%  
%  Syntax
%    alphaIndw = INDUCEDALPHA(U, wing, GAMMAw)
%    [alphaIndw, alphaIndt] = INDUCEDALPHA(U, wing, GAMMAw, tail, GAMMAt)
%
%  Input Arguments
%    U - Intensità Corrente Indisturbata [m / s]
%    wing - Matrice delle Strutture degli N x 2M Pannelli dell'Ala
%    GAMMAw - Griglia delle Circolazioni sull'Ala
%    wing - Matrice delle Strutture degli N x 2M Pannelli della Coda
%    GAMMA - Griglia delle Circolazioni sulla Coda
%
%  Output Arguments
%    alphaIndw - Vettore Contenente le Incidenze Indotte sull'Ala
%    alphaIndt - Vettore Contenente le Incidenze Indotte sulla Coda

hasTail = true;
if(nargin == 3)
    if(nargout > 1)
        error('Controllare Conformità Output con ''help inducedAlpha''')
    end
    hasTail = false;
end

% Variabili di Pannellizzazione
N = size(wing, 1);
M = size(wing, 2) / 2;

if(hasTail)
    % Variabili Pannellizzazione Coda
    Nt = size(tail, 1);
    Mt = size(tail, 2) / 2;

    % Calcolo Incidenze Indotte sull'Ala
    alphaIndw = zeros(N, 1);
    for i = 1 : 2 * M
        % Individuazione Punto di Controllo
        controlw.xc = [0; 0; 0];
        controlw.xc(1) = .25 * (.5 * (wing(end, i).x3(1) + wing(end, i).x4(1)) - .5 * (wing(1, i).x1(1) + wing(1, i).x2(1)));
        controlw.xc(3) = wing(1, i).xc(3);

        % Calcolo Velocità Indotta dai Vortici Semi-Infiniti dell'Ala
        u = 0;
        for j = 1 : N
            for k = 1 : 2 * M
                u = u + GAMMAw(j, k) * inducedVel(controlw, wing(j, k), false);
            end
        end

        % Calcolo Velocità Indotta dai Vortici Semi-Infiniti della Coda
        for j = 1 : Nt
            for k = 1 : 2 * Mt
                u = u + GAMMAt(j, k) * inducedVel(controlw, tail(j, k), false);
            end
        end

        % Calcolo Angolo d'Incidenza Indotto sull'Ala
        alphaIndw(i) = atan(u(2) / U); % Senza Diedro Normale Verticale
    end

    % Calcolo Incidenze Indotte sulla Coda
    alphaIndt = zeros(Nt, 1);
    for i = 1 : 2 * Mt
        % Individuazione Punto di Controllo
        controlt.xc = [0; 0; 0];
        controlt.xc(1) = .25 * (.5 * (tail(end, i).x3(1) + tail(end, i).x4(1)) - .5 * (tail(1, i).x1(1) + tail(1, i).x2(1)));
        controlt.xc(3) = tail(1, i).xc(3);

        % Calcolo Velocità Indotta dai Vortici Semi-Infiniti dell'Ala
        u = 0;
        for j = 1 : N
            for k = 1 : 2 * M
                u = u + GAMMAw(j, k) * inducedVel(controlt, wing(j, k), false);
            end
        end

        % Calcolo Velocità Indotta dai Vortici Semi-Infiniti della Coda
        for j = 1 : Nt
            for k = 1 : 2 * Mt
                u = u + GAMMAt(j, k) * inducedVel(controlt, tail(j, k), false);
            end
        end

        % Calcolo Angolo d'Incidenza Indotto sulla Coda
        alphaIndt(i) = atan(u(2) / U); % Senza Diedro Normale Verticale
    end
    
else
    % Calcolo Incidenze Indotte
    alphaIndw = zeros(N, 1);
    for i = 1 : 2 * M
        % Individuazione Punto di Controllo
        controlw.xc = [0; 0; 0];
        controlw.xc(1) = .25 * (.5 * (wing(end, i).x3(1) + wing(end, i).x4(1)) - .5 * (wing(1, i).x1(1) + wing(1, i).x2(1)));
        controlw.xc(3) = wing(1, i).xc(3);

        % Calcolo Velocità Indotta dai Vortici Semi-Infiniti
        u = 0;
        for j = 1 : N
            for k = 1 : 2 * M
                u = u + GAMMAw(j, k) * inducedVel(controlw, wing(j, k), false);
            end
        end

        % Calcolo Angolo d'Incidenza Indotto
        alphaIndw(i) = atan(u(2) / U); % Senza Diedro Normale Verticale
    end
end

end