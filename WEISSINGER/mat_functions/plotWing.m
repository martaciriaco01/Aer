function [x, y, z, xc, yc, zc, xt, yt, zt, xct, yct, zct] = plotWing(name, wing, tail)
%  PLOTWING - Ballielo, Ciriaco
%  Realizza il Plot 3D della Superficie dell'Ala del
%  Velivolo Assegnato. Se Fornita, Realizza anche il
%  Plot della Superficie di Coda.
%  Se Richiesti Argomenti in Output, Restituisce la
%  Griglia degli Estremi dei Pannelli delle Superfici.
%  
%  Syntax
%    PLOTWING(name, wing)
%    PLOTWING(name, wing, tail)
%    [x, y, z] = PLOTWING(name, wing)
%    [x, y, z, xt, yt, zt] = PLOTWING(name, wing, tail)
%
%  Input Arguments
%    name - Stringa Contenente il Titolo del Plot
%    wing - Matrice delle Strutture degli N x 2M Pannelli dell'Ala
%    tail - Matrice delle Strutture degli Mt x 2Mt Pannelli della Coda
%
%  Output Arguments
%    x, y, z - Matrici degli Estremi dei Pannelli dell'Ala
%    xc, yc, zc - Matrici dei Punti di Controllo dell'Ala
%    xt, yt, zt - Matrici degli Estremi dei Pannelli della Coda
%    xct, yct, zct - Matrici dei Punti di Controllo della Coda

hasTail = true;
if(nargin == 2)
    if(nargout > 6)
        error('Controllare Conformità Output con ''help plotWing''')
    end
    hasTail = false;
end

% Variabili di Pannellizzazione
N = size(wing, 1);
M = size(wing, 2) / 2;

% Griglia Estremi dei Pannelli
X = zeros(3, N + 1, 2 * M + 1);
for i = [1 : 2 : N - 1, N]
    for j = [1 : 2 : 2 * M - 1, 2 * M]
        X(:, i, j) = wing(i, j).x1;
        X(:, i, j + 1) = wing(i, j).x2;
        X(:, i + 1, j + 1) = wing(i, j).x3;
        X(:, i + 1, j) = wing(i, j).x4;
    end
end
xe = squeeze(X(1, :, :));
ye = squeeze(X(2, :, :));
ze = squeeze(X(3, :, :));

% Griglia Nodi Vortici
XV = zeros(3, N, 2 * M + 1);
for i = 1 : N
    for j = [1 : 2 : 2 * M - 1, 2 * M]
        XV(:, i, j) = wing(i, j).xv1;
        XV(:, i, j + 1) = wing(i, j).xv2;
    end
end
xv = squeeze(XV(1, :, :));
yv = squeeze(XV(2, :, :));
zv = squeeze(XV(3, :, :));

% Griglia Punti di Controllo
XC = zeros(3, N, 2 * M);
for i = 1 : N
    for j = 1 : 2 * M
        XC(:, i, j) = wing(i, j).xc;
    end
end
xC = squeeze(XC(1, :, :));
yC = squeeze(XC(2, :, :));
zC = squeeze(XC(3, :, :));

figure()
surf(ze, xe, ye)
hold on
plot3(zv, xv, yv, 'kx')
plot3(zC, xC, yC, 'ko')

if(nargout > 0)
    x = xe;
    y = ye;
    z = ze;
    xc = xC;
    yc = yC;
    zc = zC;
end

if(hasTail)
    % Variabili di Pannellizzazione
    Nt = size(tail, 1);
    Mt = size(tail, 2) / 2;

    % Griglia Estremi dei Pannelli
    Xt = zeros(3, Nt + 1, 2 * Mt + 1);
    for i = [1 : 2 : Nt - 1, Nt]
        for j = [1 : 2 : 2 * Mt - 1, 2 * Mt]
            Xt(:, i, j) = tail(i, j).x1;
            Xt(:, i, j + 1) = tail(i, j).x2;
            Xt(:, i + 1, j + 1) = tail(i, j).x3;
            Xt(:, i + 1, j) = tail(i, j).x4;
        end
    end
    xet = squeeze(Xt(1, :, :));
    yet = squeeze(Xt(2, :, :));
    zet = squeeze(Xt(3, :, :));

    % Griglia Nodi Vortici
    XVt = zeros(3, Nt, 2 * Mt + 1);
    for i = 1 : Nt
        for j = [1 : 2 : 2 * Mt - 1, 2 * Mt]
            XVt(:, i, j) = tail(i, j).xv1;
            XVt(:, i, j + 1) = tail(i, j).xv2;
        end
    end
    xvt = squeeze(XVt(1, :, :));
    yvt = squeeze(XVt(2, :, :));
    zvt = squeeze(XVt(3, :, :));

    % Griglia Punti di Controllo
    XCt = zeros(3, Nt, 2 * Mt);
    for i = 1 : Nt
        for j = 1 : 2 * Mt
            XCt(:, i, j) = tail(i, j).xc;
        end
    end
    xCt = squeeze(XCt(1, :, :));
    yCt = squeeze(XCt(2, :, :));
    zCt = squeeze(XCt(3, :, :));

    surf(zet, xet, yet)
    plot3(zvt, xvt, yvt, 'kx')
    plot3(zCt, xCt, yCt, 'ko')

    if(nargout > 6)
        xt = xet;
        yt = yet;
        zt = zet;
        xct = xCt;
        yct = yCt;
        zct = zCt;
    end
end

hold off
axis equal
view(135, 30)
title(name, 'FontSize', 16)
xlabel('z [m]', 'FontSize', 12)
ylabel('x [m]', 'FontSize', 12)
zlabel('y [m]', 'FontSize', 12)

end