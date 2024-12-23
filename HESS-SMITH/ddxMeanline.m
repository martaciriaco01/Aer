function ddx = ddxMeanline(S, xx)
%  Restituisce l'Approssimazione della Derivata Prima della Linea Media
%  (Normalizzata rispetto alla Corda) Effettuata con un Metodo di Ordine
%  di Accuratezza 2 Applicato a Nodi Equispaziati.
%  Per i Nodi Interni è Applicato uno Schema alle Differenze Finite
%  Centrate, mentre per i Nodi Estremi di Ordine Analogo.
%  
%  Syntax
%    ddx = ddxMeanline(S, xx)
%
%  Input Arguments
%    S - Vettore delle Ordinate dei Nodi della Linea Media
%    xx - Vettore delle Ascisse dei Nodi (EQUISPAZIATI)
%
%  Output Arguments
%    ddx - Valori della Derivata in xx

h = diff(xx);
if(~all(abs(h - h(1)) < 1e-12))
    error('Nodi Non Equispaziati')
end

n = length(xx);
ddx = zeros(n, 1);

% Schema al Primo Nodo (Ordine 2)
ddx(1) = (-3 * S(1) + 4 * S(2) - S(3)) / (2 * h(1));

% Schema alle Differenze Finite Centrate
for j = 2 : n - 1
    ddx(j) = (S(j + 1) - S(j - 1)) / (2 * h(1));
end

% Schema all'Ultimo Nodo (Ordine 2)
ddx(n) = (S(n - 2) - 4 * S(n - 1) + 3 * S(n)) / (2 * h(1));

end