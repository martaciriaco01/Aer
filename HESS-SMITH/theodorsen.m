%% Calcolo Angolo di Theodorsen del Profilo Alare Assegnato
% (GIII BL430 AIRFOIL)
close all
clear, clc
addpath mat_functions\

%% Importazione Dorso e Ventre del Profilo
% Necessario Disporre nel Path del File BL430.dat (Selig Format)
% Source: http://airfoiltools.com/airfoil/seligdatfile?airfoil=giiim-il
profile = importXfoilProfile('BL430.dat');

x = flipud(profile.x) - .5;
y = flipud(profile.y);
np = length(x); % Numero Punti Assegnati sul Profilo
fprintf(">> Numero Punti del Profilo Assegnati: %d\n", np)

profile = [x, y];

% Chiusura Bordo d'Uscita
profile(1, :) = [.5, 0];
profile(np, :) = [.5, 0];

% Plot Punti Assegnati
figure(1)
plot(profile(:, 1), profile(:, 2), '-o', 'LineWidth', 1.2)
hold on
title('GIII BL430 AIRFOIL', 'FontSize', 16)
xlabel('x/c', 'FontSize', 12)
ylabel('y/c', 'FontSize', 12)
grid on
axis equal
set(gca, 'XLim', [-.6, .6])

% Ipotesi (Verificata per il Profilo Assegnato):
% Stesso Numero di Punti Assegnati su Dorso e Ventre
% Di Conseguenza:
% - Se np è Pari, Bordo d'Attacco non Assegnato
% - Se np è Dispari, Bordo d'Attacco Assegnato
nm = floor(np / 2) + 1; % Numero Punti su Ciascun Side
fprintf(">> Numero Punti su Ciascun Side (Includendo LE): %d\n", nm)

% Definizione Linea del Dorso
suction = profile((nm + mod(np, 2) - 1) : np, :);
suction(1, :) = [-.5, 0]; % Aggiunta Bordo d'Attacco

% Definizione Linea del Ventre
pressure = flip(profile(1 : nm, :));
pressure(1, :) = [-.5, 0]; % Aggiunta Bordo d'Attacco

%% Individuazione Linea Media
% Ipotesi (da Verificare):
% Se, come per i Profili Naca a 4 Cifre, il Profilo è
% Definito a Partire dal Plot di una Funzione Spessore
% Ortogonalmente alla Linea Media, allora i Punti
% del Profilo Vengono Assegnati a Coppie, e ad Ogni Coppia
% Corrisponde un Punto della Linea Media Posizionato nel
% Loro Punto Medio
meanline = .5 * (suction + pressure);

% Plot Linea Media
plot(meanline(:, 1), meanline(:, 2), '-o', 'LineWidth', 1.2)
legend('Punti Profilo Alare', 'Punti Linea Media', 'FontSize', 12)
hold off

%% Verifica Validità Linea Media
% Una Possibilità è quella di Interpretare la Linea Media come
% quella Linea che Passa per Bordo d'Attacco e Bordo d'Uscita
% e che Divide il Profilo in Due Porzioni con la Stessa Area.
% Di Conseguenza, per Verificare la Sensatezza della Linea Media
% Appena Identificata, si Procede Calcolando le Aree della
% Porzione di Profilo Superiore (UP) e di quella Inferiore (DOWN)
% e Confrontandone i Valori.
% Operativamente, per Calcolare le Aree delle due Figure si
% Utilizza una Funzione che, Ricevendo come Input i Punti che
% ne Delineano il Contorno, Effettua un'Interpolazione Tramite
% Spline del Contorno Stesso e Calcola la Superficie Così Definita
% Source: https://it.mathworks.com/matlabcentral/fileexchange/69055-interpclosed/
up = [flip(suction); meanline];
down = [flip(meanline); pressure];

[len_up, area_up] = interpclosed(up(:, 1), up(:, 2));
fprintf(">> Area Porzione Superiore: %d\n", area_up)
[len_down, area_down] = interpclosed(down(:, 1), down(:, 2));
fprintf(">> Area Porzione Inferiore: %d\n", area_down)

%% Interpolazione Linea Media, Calcolo Angolo e Analisi Consistenza
% Verificato che la Linea Media è Valida, Occorre Interpolare i
% Punti che la Compongono per poi Definire la Funzione Integranda
% per il Calcolo dell'Angolo di Theodorsen.
% Si Definisce un Vettore xx Contenente un Numero n Sufficientemente
% Elevato di Nodi Dove Valutare l'Interpolazione.
% Per l'Interpolazione si Adottano 3 Diversi Metodi:
% - Spline Cubica (Not a Knot)
% - Pchip (Cubic Hermite Polynomials)
% - Interpolante Polinomiale ai Minimi Quadrati Pesati (Ordine m)
% Si Procederà con una Valutazione dell'Impatto della Scelta del
% Metodo Interpolatorio sul Risultato Finale
%
% Per la Derivazione Numerica della Linea Media, si Utilizza il Metodo
% alle Differenze Finite Centrate per i Nodi Interni e un Metodo di
% Ordine 2 agli Estremi
%
% Per Ciascun Metodo, si Calcola Numericamente l'Integrale per
% Determinare l'Angolo di Progetto Attraverso il Metodo di
% Quadratura del Trapezio.
% Essendo la Funzione Integranda Singolare agli Estremi,
% è Necessario Introdurre un Parametro epsilon, del quale si
% Studieranno gli Effetti sul Risultato, che Permetta di Distanziarsi
% dalle Singolarità nel Calcolo dell'Integrale di Theodorsen.
epsilon = [1e-5, 1e-6, 1e-7, 1e-8]; % Distanziamento dagli Estremi
n = [1e+4, 1e+5, 1e+6, 1e+7]; % Numero Nodi Interpolazione Linea Media
m = [5, 7, 9]; % Grado Fit ai Minimi Quadrati

alpha_th_spline = zeros(length(n), length(epsilon)); % Risultati Metodo Spline
alpha_th_pchip = alpha_th_spline; % Risultati Metodo Pchip
alpha_th_fit = zeros(length(m), length(n), length(epsilon)); % Risultati Metodo Minimi Quadrati

leg = {}; % Legenda Plots

for i = 1 : length(epsilon)
    fprintf("---------------- %c: 10^%d ----------------\n", char(0x03B5), log10(epsilon(i)))
    
    for j = 1 : length(n)
        fprintf("Nodi Interpolazione: 10^%d\n", log10(n(j)))
        
        xx = linspace(-.5 + epsilon(i), .5 - epsilon(i), n(j))';
        
        % Spline
        S = spline(meanline(:, 1), meanline(:, 2), xx);
        ddx = ddxMeanline(S, xx);
        f = ddx ./ sqrt(.25 - xx.^2);
        alpha_th_spline(j, i) = 180 / pi^2 * trapz(xx, f);
        fprintf(">> Spline: %c = %.4f (deg)\n", char(0x03B1), alpha_th_spline(j, i))

        % Pchip
        S = pchip(meanline(:, 1), meanline(:, 2), xx);
        ddx = ddxMeanline(S, xx);
        f = ddx ./ sqrt(.25 - xx.^2);
        alpha_th_pchip(j, i) = 180 / pi^2 * trapz(xx, f);
        fprintf(">> Pchip: %c = %.4f (deg)\n", char(0x03B1), alpha_th_pchip(j, i))

        % Minimi Quadrati Pesati
        for k = 1 : length(m)
            S = fitMeanline(meanline(:, 1) + .5, meanline(:, 2), xx + .5, m(k));
            ddx = ddxMeanline(S, xx);
            f = ddx ./ sqrt(.25 - xx.^2);
            alpha_th_fit(k, j, i) = 180 / pi^2 * trapz(xx, f);
            fprintf(">> Fit (Ord. %d): %c = %.4f (deg)\n", m(k), char(0x03B1), alpha_th_fit(k, j, i))
        end
    end

    % Valutazione Consistenza Metodo
    % Si Realizzano i Grafici dove Vengono Riportati i Valori dell'Angolo
    % di Theodorsen Ottenuti, in Relazione al Metodo Utilizzato per
    % Effettuare l'Interpolazione della Linea Media e dei Gradi di Libertà
    % Fissati, Ovvero i Parametri epsilon, n e m
    figure(2) % Metodo Spline
    semilogx(n, alpha_th_spline(:, i), '-o', 'LineWidth', 1.2)
    hold on
    figure(3) % Metodo Pchip
    semilogx(n, alpha_th_pchip(:, i), '-o', 'LineWidth', 1.2)
    hold on
    figure(4) % Metodo Minimi Quadrati
    semilogx(n, alpha_th_fit(k, :, i), '-o', 'LineWidth', 1.2)
    hold on

    leg{end + 1} = sprintf('%c = 10^{%d}', char(0x03B5), log10(epsilon(i))); % Legenda Plots
end

figure(2) % Metodo Spline
title('Interpolazione Spline', 'FontSize', 16)
legend(leg, 'FontSize', 12)
xlabel('n (Nodi Linea Media)', 'FontSize', 12)
ylabel(sprintf('%c (Angolo di Theodorsen)', char(0x03B1)), 'FontSize', 12)
grid on
hold off
figure(3) % Metodo Pchip
title('Interpolazione Pchip', 'FontSize', 16)
legend(leg, 'FontSize', 12)
xlabel('n (Nodi Linea Media)', 'FontSize', 12)
ylabel(sprintf('%c (Angolo di Theodorsen)', char(0x03B1)), 'FontSize', 12)
grid on
hold off
figure(4) % Metodo Minimi Quadrati
title(sprintf('Fit Polinomiale (Ord. %d)', m(end)), 'FontSize', 16)
legend(leg, 'FontSize', 12)
xlabel('n (Nodi Linea Media)', 'FontSize', 12)
ylabel(sprintf('%c (Angolo di Theodorsen)', char(0x03B1)), 'FontSize', 12)
grid on
hold off

%% Confronto Risultati Fit Polinomiale ai Minimi Quadrati
leg = {};
figure(5)
for k = 1 : length(m)
    semilogx(n, alpha_th_fit(k, :, i), '-o', 'LineWidth', 1.2)
    hold on
    leg{end + 1} = sprintf('m = %d', m(k));
end
title(sprintf('Confronto Ordini Fit (%c = 10^{%d})', char(0x03B5), log10(epsilon(i))), 'FontSize', 16)
legend(leg, 'FontSize', 12)
xlabel('n (Nodi Linea Media)', 'FontSize', 12)
ylabel(sprintf('%c (Angolo di Theodorsen)', char(0x03B1)), 'FontSize', 12)
grid on
hold off

