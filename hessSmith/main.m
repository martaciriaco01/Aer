%% Validazione Implementazione Metodo Hess-Smith
% Si Procede di Seguito con la Validazione dell'Implemetnazione
% del Metodo di Hess-Smith, Attraverso un Confronto dei Risultati
% Ottenuti per il Profilo NACA 0012 a Diversi Valori dell'Angolo di
% Attacco (Caso Ideale) con quelli Corrispettivi Prodotti con XFoil
close all
clear, clc
addpath mat_functions\

%% Confronto Cl
% Eseguendo lo Script Bash hessSmith.sh, si Esporta da XFoil il File
% NACA_0012.dat (200 Nodi sul Profilo), Necessario per Hess-Smith.
% Sempre con lo Script, Sfruttando la Routine "OPERi" di XFoil per
% l'Analisi Inviscida del Flusso Attorno al NACA 0012, per un
% Intervallo di Angoli d'Attacco da -4 a +8 Gradi si Ricavano i
% Risultati Utili al Successivo Confronto
%
% Per i Coefficienti Calcolati con XFoil alla Sequenza di Angoli
% d'Attacco Stabilita, si Procede con un'Ispezione del File polar.dat,
% che Contiene i Valori Tabulati di CL per ogni Alpha
results_xf = importXfoilProfile('polar.dat', 13);
alpha = results_xf.x;
cl_xf = results_xf.y;

% I Coefficienti Omologhi si Calcolano con Hess-Smith Eseguendo
% la Funzione traccia.m (Semplice Adattamento della Traccia di
% Hess-Smith Fornita per Supportare in Input una Serie di Angoli
% d'Attacco Invece di un Valore Singolo)
% !!! La Funzione è Stata Lasciata Volutamente il Più Possibile
%     Conforme alla Traccia Fornita per il Laboratorio, Infatti
%     l'Input si Limita Esclusivamente alla Serie di Alpha.
%     Sarebbe Ovviamente Possibile Generalizzarla Ulteriormente
na = length(alpha); % Numero Angoli d'Attacco
cl_hs = zeros(na, 1);
for j = 1 : na
    if(alpha(j) == 2)
        [cl_hs(j), x_hs, cp_hs] = traccia(alpha(j));
    else
        cl_hs(j) = traccia(alpha(j));
    end
end
cl_hs = round(cl_hs, 4);

% Print e Confronto Coefficienti Cl
% Si Definisce la Differenza tra i Coefficienti Ottenuti Rapportata
% al Corrispettivo Valore del Coefficiente Ricavato con XFoil
delta_cl = round(abs((cl_xf - cl_hs) ./ cl_xf), 4) * 100;
for j = 1 : na
    if(isnan(delta_cl(j)))
        delta_cl(j) = 0;
    end
end
results = table(alpha, cl_xf, cl_hs, delta_cl, 'VariableNames', {sprintf('%c (deg)', char(0x03B1)), 'Cl_xf', 'Cl_hs', sprintf('%cCl (%%)', char(0x0394))});
fprintf(">> Confronto Coefficienti Cl:\n\n")
disp(results)

figure(2)
plot(alpha, cl_xf, 'o', 'LineWidth', 1.2)
hold on
plot(alpha, cl_hs, 'LineWidth', 1.2)
title(sprintf('Curva Cl(%c) NACA 0012', char(0x03B1)), 'FontSize', 16)
legend('XFoil', 'Hess-Smith', 'Location', 'northwest', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('Cl', 'FontSize', 12)
grid on
hold off

%% Confronto Cp
% Per alpha = 2 Gradi, si Ricava il Diagramma del Coefficiente di
% Pressione Cp, Estraendo i Dati Forniti da XFoil dal File cp.dat
cp_xf = importXfoilProfile('cp.dat');

% I valori del Coefficiente di Pressione del Metodo di Hess-Smith
% sono stati Calcolati in Precedenza
% Si Realizza un Primo Plot dei Diagrammi dei Cp
figure(3)
plot(cp_xf.x(1 : 2 : end), -cp_xf.y(1 : 2 : end), 'o', 'LineWidth', 1.2)
hold on
plot(x_hs, -cp_hs, 'LineWidth', 1.2)
title(sprintf('Cp NACA 0012 (%c = 2%c)', char(0x03B1), char(0x00B0)), 'FontSize', 16)
legend('XFoil', 'Hess-Smith', 'FontSize', 12)
xlabel('x/c', 'FontSize', 12)
ylabel('-Cp', 'FontSize', 12)
grid on
hold off

% Si Vuole Approfondire il Confronto tra i Diagrammi.
% Perciò, si Procede Definendo una Differenza Relativa
% che Faciliti il Confronto tra i Due
%
% Ipotesi (Verificata per il Profilo NACA 0012):
% Stesso Numero di Punti Assegnati su Dorso e Ventre
% Di Conseguenza:
% - Se np è Pari, Bordo d'Attacco non Assegnato
% - Se np è Dispari, Bordo d'Attacco Assegnato
np = length(cp_xf.x); % XFoil Assegna il Cp sui Nodi del Profilo
nm = ceil(np / 2); % Numero Nodi su Ciascun Side
x_xf_suction = flipud(cp_xf.x(1 : nm)); % Ascisse Nodi Dorso
cp_xf_suction = flipud(cp_xf.y(1 : nm));
x_xf_pressure = cp_xf.x(nm - mod(np, 2) + 1 : np); % Ascisse Nodi Ventre
cp_xf_pressure = cp_xf.y(nm - mod(np, 2) + 1 : np);

np = length(x_hs); % Hess-Smith Calcola Cp nei Punti Medi dei Pannelli
nm = ceil(np / 2); % Numero Punti Medi su Ciascun Side
x_hs_suction = x_hs(nm - mod(np, 2) + 1 : np); % Ascisse Punti Medi Dorso
cp_hs_suction = cp_hs(nm - mod(np, 2) + 1 : np);
x_hs_pressure = flipud(x_hs(1 : nm)); % Ascisse Punti Medi Ventre
cp_hs_pressure = flipud(cp_hs(1 : nm));

% Dal Momento che non c'è Diretta Corrispondenza tra le Ascisse
% Rispetto alle quali Vengono Assegnati i Valori di Cp dai due
% Strumenti, per Rendere Confrontabili i Diagrammi si Effettua
% un'Interpolazione Lineare dei Dati
x_suction = linspace(max(x_xf_suction(1), x_hs_suction(1)), min(x_xf_suction(end), x_hs_suction(end)), 1e+3)';
cp_xf_suction = interp1(x_xf_suction, cp_xf_suction, x_suction);
cp_hs_suction = interp1(x_hs_suction, cp_hs_suction, x_suction);
x_pressure = linspace(max(x_xf_pressure(1), x_hs_pressure(1)), min(x_xf_pressure(end), x_hs_pressure(end)), 1e+3)';
cp_xf_pressure = interp1(x_xf_pressure, cp_xf_pressure, x_pressure);
cp_hs_pressure = interp1(x_hs_pressure, cp_hs_pressure, x_pressure);

% Calcolo Differenza Relativa Cp
% Si Definisce la Differenza dei Diagrammi del Cp di Dorso e Ventre
% (è già una Grandezza Relativa, Rapportata al Cp(LE) = 1)
delta_cp_suction = abs(cp_xf_suction - cp_hs_suction);
delta_cp_pressure = abs(cp_xf_pressure - cp_hs_pressure);

% Plot Differenze
figure(4)
plot(x_suction, delta_cp_suction, 'LineWidth', 1.2)
hold on
plot(x_pressure, delta_cp_pressure, 'LineWidth', 1.2)
axis equal
title(sprintf('Differenza Relativa %cCp(%%)', char(0x0394)), 'FontSize', 16)
legend('Dorso', 'Ventre', 'FontSize', 12)
xlabel('x/c', 'FontSize', 12)
ylabel(sprintf('%cCp(%%)', char(0x0394)), 'FontSize', 12)
grid on
hold off