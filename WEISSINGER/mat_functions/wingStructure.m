function [wing, S] = wingStructure(c, b, taper1, sweep1, M, N, b2, taper2, sweep2, off)
%  WINGSTRUCTURE - Ballielo, Ciriaco
%  Restituisce la Pannellizzazione dell'Ala del Velivolo
%  a Partire dai Parametri Geometrici in input.
%  Supporta Velivoli la cui Semiala Presenta 2 Diversi
%  Tratti con Differente Sweep Angle e Taper Ratio.
%  Inoltre, Calcola la Superficie dell'Ala.
%  
%  Syntax
%    [wing, S] = WINGSTRUCTURE(c, b, taper1, taper2, sweep, M, N, b2, taper2, sweep2, off)
%
%  Input Arguments
%    c - Corda Aerodinamica alla Radice [m]
%    b - Apertura Alare Completa [m]
%    taper1 - Taper Ratio del Primo Tratto (c(b/2 - b2) / c(0))
%    sweep1 - Angolo di Freccia del Primo Tratto (A Metà Corda) [deg]
%    M - Numero Settori Longitudinali Discretizzazione (Semi-Ala)
%    N - Numero Pannelli Discretizzazione (PARI)
%    b2 - Semi Estensione del Secondo Tratto [m]
%    taper2 - Taper Ratio del Secondo Tratto (c(b/2) / c(b/2 - b2))
%    sweep2 - Angolo di Freccia del Secondo Tratto (A Metà Corda) [deg]
%    off - Vettore Offset Ala [m]
%
%  Output Arguments
%    wing - Matrice Contenente le Strutture degli N x 2M Pannelli
%           Ciascuna Struttura Contiene:
%           x1, x2, x3, x4 - Vettori Posizione degli Estremi
%           xv1, xv2 - Vettori Posizione dei Nodi dei Vortici
%           xc - Vettore Posizione del Punto di Controllo
%   S - Superficie Ala [m^2]

is2steps = true;

if(nargin == 6)
    is2steps = false;
    off = zeros(3, 1); % Nessun Offset
elseif(nargin == 7)
    is2steps = false;
    off = b2;
elseif(nargin == 9)
    off = zeros(3, 1);
elseif(nargin ~= 10)
    error('Controllare Conformità Input con ''help wingStructure''')
end

% Verifica N Pari
if(mod(N, 2) ~= 0)
    error('N Dispari')
end

sweep1 = sweep1 / 180 * pi;

if(is2steps)
    sweep2 = sweep2 / 180 * pi;

    if(b2 < b / (4 * M) || b2 >= b * (2 * M - 1) / (4 * M))
        error('Valore b2 Non Idoneo')
    end

    % Divisione Tratti
    [~, m] = min(abs(linspace(0, b / 2, M + 1) - b / 2 + b2));
    m = m - 1;

    % Griglia Estremi dei Pannelli
    % Totale Estremi: (2M + 1)(N + 1)
    % Griglia Nodi Vortici
    % Totale Nodi Vortici: (2M + 1)(N + 1)
    % Griglia Punti di Controllo
    % Totale Punti di Controllo: 2M * N
    d = b2 / (M - m);
    B = b / 2 - d / 2;
    z = [b / 2 : -d : b / 2 - b2 + d, linspace(b / 2 - b2, b2 - b / 2, 2 * m + 1), -b / 2 + b2 - d : -d : -b / 2];
    D = (b / 2 - b2) / (2 * m);
    zc = [B : -d : B - b2 + d, linspace(b / 2 - b2 - D, -b / 2 + b2 + D, 2 * m), -B + b2 - d : -d : -B];
    x = zeros(N + 1, 2 * M + 1);
    xv = zeros(N, 2 * M + 1);
    xc = zeros(N, 2 * M);
    c2 = taper1 * c; % Corda al Cambio di Tratto
    for j = M - m + 1 : M + m + 1
        t1 = 1 - ((1 - taper1) / (b / 2 - b2)) * abs(z(j)); % Contributo Taper Ratio 1
        s1 = abs(z(j)) * tan(sweep1); % Contributo Sweep Angle 1
        x(:, j) = c / 2 + s1 + t1 * linspace(-c / 2, c / 2, N + 1)';
        d = c / (4 * N);
        xv(:, j) = c / 2 + s1 + t1 * linspace(-c / 2 + d, c / 2 - 3 * d, N)';
        t = 1 - ((1 - taper1) / (b / 2 - b2)) * abs(zc(j)); % Contributo Taper Ratio 1
        s = abs(zc(j)) * tan(sweep1); % Contributo Sweep Angle 1
        if(j == M + m + 1)
            break
        end
        xc(:, j) = c / 2 + s + t * linspace(-c / 2 + 3 * d, c / 2 - d, N)';
    end
    for j = 1 : M - m
        t2 = 1 - ((1 - taper2) / b2) * abs(z(j) - b / 2 + b2); % Contributo Taper Ratio 2
        s2 = abs(z(j) - b / 2 + b2) * tan(sweep2); % Contributo Sweep Angle 2
        x(:, j) = c / 2 + s1 + s2 + t2 * linspace(-c2 / 2, c2 / 2, N + 1)';
        d = c2 / (4 * N);
        xv(:, j) = c / 2 + s1 + s2 + t2 * linspace(-c2 / 2 + d, c2 / 2 - 3 * d, N)';
        t2 = 1 - ((1 - taper2) / b2) * abs(zc(j) - b / 2 + b2); % Contributo Taper Ratio 2
        s2 = abs(zc(j) - b / 2 + b2) * tan(sweep2); % Contributo Sweep Angle 2
        xc(:, j) = c / 2 + s1 + s2 + t2 * linspace(-c2 / 2 + 3 * d, c2 / 2 - d, N)';
    end
    x(:, M + m + 2 : 2 * M + 1) = flip(x(:, 1 : M - m), 2);
    xv(:, M + m + 2 : 2 * M + 1) = flip(xv(:, 1 : M - m), 2);
    xc(:, M + m + 1 : 2 * M) = flip(xc(:, 1 : M - m), 2);

    c3 = taper2 * c2; % Corda all'Estremità dell'Ala
    S = c * (b / 2 - b2) + c2 * b / 2 + c3 * b2;
else
    % Griglia Estremi dei Pannelli
    % Totale Estremi: (2M + 1) * (N + 1)
    % Griglia Nodi Vortici
    % Totale Nodi Vortici: (2M + 1) * N
    % Griglia Punti di Controllo
    % Totale Punti di Controllo: 2M * N
    B = b * (2 * M - 1) / (4 * M);
    z = linspace(b / 2, -b / 2, 2 * M + 1);
    zc = linspace(B, -B, 2 * M);
    x = zeros(N + 1, 2 * M + 1);
    xv = zeros(N, 2 * M + 1);
    xc = zeros(N, 2 * M);
    for j = 1 : 2 * M + 1
        t = 1 - ((1 - taper1) / (b / 2)) * abs(z(j)); % Contributo Taper Ratio
        s = abs(z(j)) * tan(sweep1); % Contributo Sweep Angle
        x(:, j) = c / 2 + s + t * linspace(-c / 2, c / 2, N + 1)';
        d = c / (4 * N);
        xv(:, j) = c / 2 + s + t * linspace(-c / 2 + d, c / 2 - 3 * d, N)';
        if(j == 2 * M  + 1)
            break
        end
        t = 1 - ((1 - taper1) / (b / 2)) * abs(zc(j)); % Contributo Taper Ratio
        s = abs(zc(j)) * tan(sweep1); % Contributo Sweep Angle
        xc(:, j) = c / 2 + s + t * linspace(-c / 2 + 3 * d, c / 2 - d, N)';
    end

    c2 = taper1 * c; % Corda all'Estremità dell'Ala
    S = (c + c2) * b / 2;
end

% Aggiunta Offset
x = x + off(1);
xv = xv + off(1);
xc = xc + off(1);
y = off(2);
z = z + off(3);
zc = zc + off(3);

% Struttura Pannello
d = zeros(3, 1);
panel = struct('x1', d, 'x2', d, 'x3', d, 'x4', d,'xv1', d, 'xv2', d, 'xc', d);

% Struttura Ala (Griglia Pannelli)
wing = repmat(panel, N, 2 * M);

for i = 1 : N
    for j = 1 : 2 * M
        wing(i, j).x1 = [x(i, j); y; z(j)];
        wing(i, j).x2 = [x(i, j + 1); y; z(j + 1)];
        wing(i, j).x3 = [x(i + 1, j + 1); y; z(j + 1)];
        wing(i, j).x4 = [x(i + 1, j); y; z(j)];
        wing(i, j).xv1 = [xv(i, j); y; z(j)];
        wing(i, j).xv2 = [xv(i, j + 1); y; z(j + 1)];
        wing(i, j).xc = [xc(i, j); y; zc(j)];
    end
end

end