function [wing, tail, GAMMAw, GAMMAt] = computeGamma(alpha, U, wing, tail)
%  COMPUTEGAMMA - Ballielo, Ciriaco
%  Restituisce la Struttura dell'Ala del Velivolo
%  con l'Aggiunta del Campo Contenente il Valore della
%  Gamma (Circolazione) Associata ad Ogni Singolo Pannello.
%  Supporta anche l'Interazione col Piano di Coda.
%  Se Richiesto, Restituisce anche la Griglia delle Gamma
%  
%  Syntax
%    wing = COMPUTEGAMMA(alpha, U, wing)    
%    [wing, GAMMAw] = COMPUTEGAMMA(alpha, U, wing)
%    [wing, tail] = COMPUTEGAMMA(alpha, U, wing, tail)
%    [wing, tail, GAMMAw, GAMMAt] = COMPUTEGAMMA(alpha, U, wing, tail)
%
%  Input Arguments
%    wing - Matrice delle Strutture degli N x 2M Pannelli
%    tail - Matrice delle Strutture degli Nt x 2Mt Pannelli
%    alpha - Angolo d'Attacco [deg]
%    U - Intensità Corrente Indisturbata [m / s]
%
%  Output Arguments
%    wing - Matrice delle Strutture degli N x 2M Pannelli
%           Ciascuna Struttura Aggiornata Contiene:
%           x1, x2, x3, x4 - Vettori Posizione degli Estremi
%           xv1, xv2 - Vettori Posizione dei Nodi dei Vortici
%           xc - Vettore Posizione del Punto di Controllo
%           gamma - Intensità del Vortice a Ferro di Cavallo
%    tail - Matrice delle Strutture degli Mt x 2Mt Pannelli
%    GAMMAw - Griglia Contenente le Circolazioni dell'Ala
%    GAMMAt - Griglia Contenente le Circolazioni della Coda

hasTail = true;
if(nargin == 3)
    if(nargout > 2)
        error('Controllare Conformità Output con ''help computeGamma''')
    end
    hasTail = false;
end

% Conversione a Radianti
alpha = alpha / 180 * pi;

% Corrente Indisturbata
U = U * [cos(alpha); sin(alpha)];

% Variabili di Pannellizzazione
N = size(wing, 1);
M = size(wing, 2) / 2;

% Reshape Matrice Struttura Pannelli a Vettore Colonna
wingPanels = N * 2 * M;
wing = reshape(wing', wingPanels, 1);

if(hasTail)
    % Variabili Pannellizzazione Coda
    Nt = size(tail, 1);
    Mt = size(tail, 2) / 2;

    % Reshape Matrice Struttura Pannelli a Vettore Colonna
    tailPanels = Nt * 2 * Mt;
    tail = reshape(tail', tailPanels, 1);

    % Impostazione Sistema A * gamma = b
    % Costruzione Matrice dei Coefficienti A
    % (Va Considerata l'Induzione Reciproca dei Piani)
    A = zeros(wingPanels + tailPanels);
    for i = 1 : wingPanels
        for j = 1 : wingPanels % Sottomatrice Aww
            u = inducedVel(wing(i), wing(j)); % Indotta da j su i
            A(i, j) = u(2);
        end
        for j = 1 : tailPanels % Sottomatrice Atw
            u = inducedVel(wing(i), tail(j));
            A(i, wingPanels + j) = u(2);
        end
    end
    for i = 1 : tailPanels
        for j = 1 : wingPanels % Sottomatrice Awt
            u = inducedVel(tail(i), wing(j));
            A(wingPanels + i, j) = u(2);
        end
        for j = 1 : tailPanels
            u = inducedVel(tail(i), tail(j));
            A(wingPanels + i, wingPanels + j) = u(2);
        end
    end
    % Costruzione del Vettore dei Termini Noti b
    b = -U(2) * ones(wingPanels + tailPanels, 1); % Senza Diedro Normale Verticale

    % Risoluzione Sistema Lineare
    gammaVector = linsolve(A, b);

    % Assegnazione Circolazione al Singolo Pannello dell'Ala
    for i = 1 : wingPanels
        wing(i).gamma = gammaVector(i);
    end
    wing = reshape(wing, 2 * M, N)';
    GAMMAw = reshape(gammaVector(1 : wingPanels), 2 * M, N)';

    % Assegnazione Circolazione al Singolo Pannello della Coda
    for i = 1 : tailPanels
        tail(i).gamma = gammaVector(wingPanels + i);
    end
    tail = reshape(tail, 2 * Mt, Nt)';
    GAMMAt = reshape(gammaVector(wingPanels + 1 : wingPanels + tailPanels), 2 * Mt, Nt)';

else
    % Impostazione Sistema A * gamma = b
    % Costruzione Matrice dei Coefficienti A
    A = zeros(wingPanels);
    for i = 1 : wingPanels
        for j = 1 : wingPanels
            u = inducedVel(wing(i), wing(j)); % Indotta da j su i
            A(i, j) = u(2); % Senza Diedro Normale Verticale
        end
    end

    % Costruzione del Vettore dei Termini Noti b
    b = -U(2) * ones(wingPanels, 1); % Senza Diedro Normale Verticale

    % Risoluzione Sistema Lineare
    gammaVector = linsolve(A, b);

    % Assegnazione Circolazione al Singolo Pannello
    for i = 1 : wingPanels
        wing(i).gamma = gammaVector(i);
    end

    wing = reshape(wing, 2 * M, N)';
    GAMMAw = reshape(gammaVector, 2 * M, N)';
    % Output GAMMAw se Richiesto
    tail = GAMMAw;
end

end