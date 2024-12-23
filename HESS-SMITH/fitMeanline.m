function p = fitMeanline(x, y, xx, m)
%  Restituisce l'Approssimazione della Linea Media di un Profilo Alare
%  (Normalizzata rispetto alla Corda) Effettuata con un Metodo Polinomiale
%  ai Minimi Quadrati Pesato.
%  Come Funzione Peso si Utilizza, nell'Intervallo [0, 1], la Parabola
%  w(x) = (2x - 1)^6, Unitaria in 0 e 1 e Nulla al Centro dell'Intervallo.
%  Come Condizioni al Contorno del Polinomio vengono Imposti
%  il Passaggio per LE e quello per TE.
%  
%  Syntax
%    S = fitMeanline(x, y, xx, n)
%
%  Input Arguments
%    x - Vettore delle Ascisse dei Nodi della Linea Media
%    y - Vettore delle Ordinate dei Nodi della Linea Media
%    xx - Vettore dei Punti in cui Valutare l'Approssimazione
%    m - Grado del Polinomio Interpolante
%
%  Output Arguments
%    p - Valori del Polinomio in xx

n = length(x) - 1;
if(n < m)
    error('Grado Interpolante Maggiore o Uguale al Numero di Nodi')
end
if(size(x, 1) < size(x, 2))
    x = x';
end
if(size(y, 1) < size(y, 2))
    y = y';
end
if(size(xx, 1) < size(xx, 2))
    xx = xx';
end

% Si Definisce la Funzione Peso (Più Rilevante in LE e TE)
w = @(x) (2 * x - 1).^6;

% Si Definisce il Vettore dei Coefficienti Incogniti del Polinomio
% Interpolante pm(x) = a_0 + a_1 * x + a_2 * x^2 + ... + a_m * x^m
a = zeros(m + 1, 1);

% Si Imposta il Sistema Lineare A * a = q
% Condizione al Bordo d'Attacco: pm(0) = 0 => a_0 = 0
% Condizioni per i = 1, ..., n => Teoria dei Minimi Quadrati
V = zeros(n + 1, m + 1); % Matrice di Vandermonde
for i = 1 : n
    for j = 1 : m + 1
        V(i, j) = x(i)^(j - 1);
    end
end
A = V' * (w(x) .* V);
A = A(2 : m + 1, 2 : m + 1);
q = V' * (w(x) .* y);
q = q(2 : m + 1);

% Condizione al Bordo d'Uscita: pm(1) = 0 => a_1 + ... + a_m = 0
A(m, :) = ones(1, m);
q(m) = 0;

% Risoluzione Sistema Lineare
a(2 : m + 1) = linsolve(A, q);

% Definizione Polinomio
pm = @(x) x.^(0 : m) * a;

% Valutazione Polinomio in xx
p = pm(xx);

end
