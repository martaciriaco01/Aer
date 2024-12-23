%% Metodo di Weissinger
% Si Utilizza l'Implementazione del Metodo di Weissinger per
% Studiare le Caratteristiche Aerodinamiche dell'Ala del
% Cessna 172 Skyhawk e del Tecnam P2008
close all
clear, clc
addpath mat_functions\

%% Parametrizzazione e Pannellizzazione Cessna 172 Skyhawk
% Parametrizzazione Ala
c_C172 = 1.65; % [m] Corda alla Radice
b_C172 = 11; % [m] Apertura Alare Completa
taper1_C172 = 1; % Taper Ratio Primo Tratto
sweep1_C172 = 0; % [deg] Sweep Angle Primo Tratto (Centro Corda)
b2_C172 = 3; % [m] Estensione Secondo Tratto
taper2_C172 = 2 / 3; % Taper Ratio Secondo Tratto
sweep2_C172 = -2; % [deg] Sweep Angle Secondo Tratto (Centro Corda)
M_C172 = 10; % Settori Longitudinali Semiala
N_C172 = 10; % Settori Lungo la Corda

% Creazione Struttura Pannelli e Calcolo Superficie
[wing_C172, S_C172] = wingStructure(c_C172, b_C172, taper1_C172, sweep1_C172, M_C172, N_C172, b2_C172, taper2_C172, sweep2_C172);
%[wing_C172, S_C172] = wingStructure(c_C172, b_C172, taper1_C172, sweep1_C172, M_C172, N_C172);
AR_C172 = b_C172^2 / S_C172;
fprintf(">> Aspect Ratio Ala Cessna 172: %.4f\n", AR_C172)

% Plot 3D
[~, ~, z_C172, ~, ~, zc_C172] = plotWing('Cessna 172 Skyhawk', wing_C172);
xlim([-6, 6])
ylim([-.5, 2])
zlim([-.5, .5])

%% Parametrizzazione e Pannellizzazione Tecnam P2008
c_P2008 = 1.45; % [m] Corda alla Radice
b_P2008 = 9; % [m] Apertura Alare Completa
taper1_P2008 = 1; % Taper Ratio Primo Tratto
sweep1_P2008 = 0; % [deg] Sweep Angle Primo Tratto (Centro Corda)
b2_P2008 = 1.7; % [m] Estensione Secondo Tratto
taper2_P2008 = 3 / 4; % Taper Ratio Secondo Tratto
sweep2_P2008 = -6; % [deg] Sweep Angle Secondo Tratto (Centro Corda)
M_P2008 = 10; % Settori Longitudinali Semiala
N_P2008 = 10; % Settori Lungo la Corda

% Creazione Struttura Pannelli e Calcolo Superficie
[wing_P2008, S_P2008] = wingStructure(c_P2008, b_P2008, taper1_P2008, sweep1_P2008, M_P2008, N_P2008, b2_P2008, taper2_P2008, sweep2_P2008);
%[wing_P2008, S_P2008] = wingStructure(c_P2008, b_P2008, taper1_P2008, sweep1_P2008, M_P2008, N_P2008);
AR_P2008 = b_P2008^2 / S_P2008;
fprintf(">> Aspect Ratio Ala Tecnam P2008: %.4f\n", AR_P2008)

% Plot 3D
[~, ~, z_P2008, ~, ~, zc_P2008] = plotWing('Tecnam P2008', wing_P2008);
xlim([-5, 5])
ylim([-.5, 2])
zlim([-.5, .5])

%% Parametri Aerodinamici
alpha = -6 : 8; % [deg] Sequenza Angoli d'Attacco
lenAlpha = length(alpha);
U = 1; % [m/s] Intensità Corrente
rho = 1.225; % [kg/m^3] Densità Aria

%% Calcolo Caratteristiche Aerodinamiche Ala Cessna 172 Skyhawk
L_C172 = zeros(lenAlpha, 1);
CL_C172 = zeros(lenAlpha, 1);
Di_C172 = zeros(lenAlpha, 1);
CDi_C172 = zeros(lenAlpha, 1);
q_C172 = .5 * rho * U^2 * S_C172;
for i = 1 : lenAlpha
    % Determinazione Circolazione Gamma
    [wing_C172, GAMMAw_C172] = computeGamma(alpha(i), U, wing_C172);
    if(alpha(i) == 2)
        GAMMAw2_C172 = GAMMAw_C172(1, :);
    end
    
    % Calcolo Portanza
    L_C172(i) = rho * U * sum(sum(GAMMAw_C172) * abs(diff(z_C172(1, :)))');

    % Calcolo CL
    CL_C172(i) = L_C172(i) / q_C172;

    % Calcolo Resistenza Indotta
    alphaInd = inducedAlpha(U, wing_C172, GAMMAw_C172);
    Di_C172(i) = -rho * U * sum(sum(GAMMAw_C172) * (abs(diff(z_C172(1, :)))' .* sin(alphaInd))); % Si Prende Positiva
    
    % Calcolo CDi
    CDi_C172(i) = Di_C172(i) / q_C172;
end

%% Calcolo Caratteristiche Aerodinamiche Ala Tecnam P2008
L_P2008 = zeros(lenAlpha, 1);
CL_P2008 = zeros(lenAlpha, 1);
Di_P2008 = zeros(lenAlpha, 1);
CDi_P2008 = zeros(lenAlpha, 1);
q_P2008 = .5 * rho * U^2 * S_P2008;
for i = 1 : lenAlpha
    % Determinazione Circolazione Gamma
    [wing_P2008, GAMMAw_P2008] = computeGamma(alpha(i), U, wing_P2008);
    if(alpha(i) == 2)
        GAMMAw2_P2008 = GAMMAw_P2008(1, :);
    end
    
    % Calcolo Portanza
    L_P2008(i) = rho * U * sum(sum(GAMMAw_P2008) * abs(diff(z_P2008(1, :)))');

    % Calcolo CL
    CL_P2008(i) = L_P2008(i) / q_P2008;

    % Calcolo Resistenza Indotta
    alphaInd = inducedAlpha(U, wing_P2008, GAMMAw_P2008);
    Di_P2008(i) = -rho * U * sum(sum(GAMMAw_P2008) * (abs(diff(z_P2008(1, :)))' .* sin(alphaInd))); % Si Prende Positiva
    
    % Calcolo CDi
    CDi_P2008(i) = Di_P2008(i) / q_P2008;
end

%% Confronto (Solo Ali)
% Plot CL-alpha
figure(3)
plot(alpha, CL_C172, '-o', 'LineWidth', 1.2)
hold on
plot(alpha, CL_P2008, '-o', 'LineWidth', 1.2)
title(sprintf('Curva CL(%c) Solo Ala', char(0x03B1)), 'FontSize', 16)
legend('Cessna 172', 'Tecnam P2008', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CL', 'FontSize', 12)
grid on
hold off

% Plot CDi-alpha
figure(4)
plot(alpha, CDi_C172, '-o', 'LineWidth', 1.2)
hold on
plot(alpha, CDi_P2008, '-o', 'LineWidth', 1.2)
title(sprintf('Curva CDi(%c) Solo Ala', char(0x03B1)), 'FontSize', 16)
legend('Cessna 172', 'Tecnam P2008', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CDi', 'FontSize', 12)
grid on
hold off

% Plot Polare
figure(5)
plot(CDi_C172, CL_C172, '-o', 'LineWidth', 1.2)
hold on
plot(CDi_P2008, CL_P2008, '-o', 'LineWidth', 1.2)
title('Polare Solo Ala', 'FontSize', 16)
legend('Cessna 172', 'Tecnam P2008', 'FontSize', 12)
xlabel('CDi', 'FontSize', 12)
ylabel('CL', 'FontSize', 12)
grid on
hold off

%% Confronto Distribuzione Ellittica (Primo Settore Longitudinale N = 1, alpha = 2 deg)
maxGamma_C172 = max(GAMMAw2_C172);
ellipticGamma_C172 = maxGamma_C172 * sqrt(1 - (2 / b_C172 * zc_C172(1, :)).^2);

maxGamma_P2008 = max(GAMMAw2_P2008);
ellipticGamma_P2008 = maxGamma_P2008 * sqrt(1 - (2 / b_P2008 * zc_P2008(1, :)).^2);

figure(6)
plot(zc_C172(1, :), GAMMAw2_C172, '-o', 'LineWidth', 1.2)
hold on
plot(zc_P2008(1, :), GAMMAw2_P2008, '-o', 'LineWidth', 1.2)
plot(zc_C172(1, :), ellipticGamma_C172, '--', 'Color', '#808080', 'LineWidth', 1.2)
plot(zc_P2008(1, :), ellipticGamma_P2008, '--', 'Color', '#808080', 'LineWidth', 1.2)
title(sprintf('Dist. Longitudinale %c (%c = 2%c, N = 1)', char(0x0393), char(0x03B1), char(0x00B0)), 'FontSize', 16)
legend('Cessna 172', 'Tecnam P2008', sprintf('%c Ellittica', char(0x0393)), 'FontSize', 12)
xlabel('z [m]', 'FontSize', 12)
ylabel(sprintf('%c [m^2/s]', char(0x0393)))
grid on
hold off

%% Introduzione Effetto dei Piani di Coda
% Si Implementa l'Interazione dell'Ala dei due Velivoli
% Con i Rispettivi Piani di Coda
close all

%% Parametrizzazione e Pannellizzazione Coda Cessna 172 Skyhawk
% Parametrizzazione Coda
ct_C172 = 1.35; % [m] Corda alla Radice
bt_C172 = 3.4; % [m] Apertura Alare Completa
tapert_C172 = 5 / 9; % Taper Ratio
sweept_C172 = 0; % [deg] Sweep Angle (Centro Corda)
off_C172 = [4.3, -.7, 0]; % [m] Offset Coda
Mt_C172 = 6; % Settori Longitudinali Semiala
Nt_C172 = 6; % Settori Lungo la Corda

% Creazione Struttura Pannelli e Calcolo Superficie
[tail_C172, St_C172] = wingStructure(ct_C172, bt_C172, tapert_C172, sweept_C172, Mt_C172, Nt_C172, off_C172);
ARt_C172 = bt_C172^2 / St_C172;
fprintf(">> Aspect Ratio Coda Cessna 172: %.4f\n", ARt_C172)

% Plot 3D
[~, ~, z_C172, ~, ~, zc_C172, ~, ~, zt_C172, ~, ~, zct_C172] = plotWing('Cessna 172 Skyhawk', wing_C172, tail_C172);
xlim([-6, 6])
ylim([-.5, 6])
zlim([-1, .5])

%% Parametrizzazione e Pannellizzazione Coda Tecnam P2008
% Parametrizzazione Coda
ct_P2008 = .7; % [m] Corda alla Radice
bt_P2008 = 3; % [m] Apertura Alare Completa
tapert_P2008 = 1; % Taper Ratio
sweept_P2008 = 0; % [deg] Sweep Angle (Centro Corda)
off_P2008 = [4, -.7, 0]; % [m] Offset Coda
Mt_P2008 = 6; % Settori Longitudinali Semiala
Nt_P2008 = 6; % Settori Lungo la Corda

% Creazione Struttura Pannelli e Calcolo Superficie
[tail_P2008, St_P2008] = wingStructure(ct_P2008, bt_P2008, tapert_P2008, sweept_P2008, Mt_P2008, Nt_P2008, off_P2008);
ARt_P2008 = bt_P2008^2 / St_P2008;
fprintf(">> Aspect Ratio Coda Tecnam P2008: %.4f\n", ARt_P2008)

% Plot 3
[~, ~, z_P2008, ~, ~, zc_P2008, ~, ~, zt_P2008, ~, ~, zct_P2008] = plotWing('Tecnam P2008', wing_P2008, tail_P2008);
xlim([-5, 5])
ylim([-.5, 5])
zlim([-1, .5])

%% Effetto Piani di Coda Cessna 172
L_C172_tail = zeros(lenAlpha, 1);
CL_C172_tail = zeros(lenAlpha, 1);
Di_C172_tail = zeros(lenAlpha, 1);
CDi_C172_tail = zeros(lenAlpha, 1);
q_C172 = .5 * rho * U^2 * S_C172;
for i = 1 : lenAlpha
    % Determinazione Circolazione Gamma
    [wing_C172, tail_C172, GAMMAw_C172, GAMMAt_C172] = computeGamma(alpha(i), U, wing_C172, tail_C172);
    
    % Calcolo Portanza Ala
    Lw_C172 = rho * U * sum(sum(GAMMAw_C172) * abs(diff(z_C172(1, :)))');

    % Calcolo Portanza Coda
    Lt_C172 = rho * U * sum(sum(GAMMAt_C172) * abs(diff(zt_C172(1, :)))');

    % Calcolo Portanza Totale
    L_C172_tail(i) = Lw_C172 + Lt_C172;

    % Calcolo CL
    CL_C172_tail(i) = L_C172_tail(i) / q_C172; % Comunque Riferito alla Superficie dell'Ala

    % Calcolo Resistenza Indotta Ala
    [alphaIndw, alphaIndt] = inducedAlpha(U, wing_C172, GAMMAw_C172, tail_C172, GAMMAt_C172);
    Diw_C172 = -rho * U * sum(sum(GAMMAw_C172) * (abs(diff(z_C172(1, :)))' .* sin(alphaIndw))); % Si Prende Positiva

    % Calcolo Resistenza Indotta Coda
    Dit_C172 = -rho * U * sum(sum(GAMMAt_C172) * (abs(diff(zt_C172(1, :)))' .* sin(alphaIndt))); % Si Prende Positiva
    
    % Calcolo Resistenza Indotta Totale
    Di_C172_tail(i) = Diw_C172 + Dit_C172;

    % Calcolo CDi
    CDi_C172_tail(i) = Di_C172_tail(i) / q_C172; % Comunque Riferito alla Superficie dell'Ala
end

%% Effetto Piani di Coda Tecnam P2008
L_P2008_tail = zeros(lenAlpha, 1);
CL_P2008_tail = zeros(lenAlpha, 1);
Di_P2008_tail = zeros(lenAlpha, 1);
CDi_P2008_tail = zeros(lenAlpha, 1);
q_P2008 = .5 * rho * U^2 * S_P2008;
for i = 1 : lenAlpha
    % Determinazione Circolazione Gamma
    [wing_P2008, tail_P2008, GAMMAw_P2008, GAMMAt_P2008] = computeGamma(alpha(i), U, wing_P2008, tail_P2008);
    
    % Calcolo Portanza Ala
    Lw_P2008 = rho * U * sum(sum(GAMMAw_P2008) * abs(diff(z_P2008(1, :)))');

    % Calcolo Portanza Coda
    Lt_P2008 = rho * U * sum(sum(GAMMAt_P2008) * abs(diff(zt_P2008(1, :)))');

    % Calcolo Portanza Totale
    L_P2008_tail(i) = Lw_P2008 + Lt_P2008;

    % Calcolo CL
    CL_P2008_tail(i) = L_P2008_tail(i) / q_P2008; % Comunque Riferito alla Superficie dell'Ala

    % Calcolo Resistenza Indotta Ala
    [alphaIndw, alphaIndt] = inducedAlpha(U, wing_P2008, GAMMAw_P2008, tail_P2008, GAMMAt_P2008);
    Diw_P2008 = -rho * U * sum(sum(GAMMAw_P2008) * (abs(diff(z_P2008(1, :)))' .* sin(alphaIndw))); % Si Prende Positiva

    % Calcolo Resistenza Indotta Coda
    Dit_P2008 = -rho * U * sum(sum(GAMMAt_P2008) * (abs(diff(zt_P2008(1, :)))' .* sin(alphaIndt))); % Si Prende Positiva
    
    % Calcolo Resistenza Indotta Totale
    Di_P2008_tail(i) = Diw_P2008 + Dit_P2008;

    % Calcolo CDi
    CDi_P2008_tail(i) = Di_P2008_tail(i) / q_P2008; % Comunque Riferito alla Superficie dell'Ala
end

%% Confronto Cessna 172 (Con e Senza Piani di Coda)
% Plot CL-alpha
figure(3)
plot(alpha, CL_C172, '-o', 'LineWidth', 1.2)
hold on
plot(alpha, CL_C172_tail, '-o', 'LineWidth', 1.2)
title(sprintf('Curva CL(%c) Cessna 172', char(0x03B1)), 'FontSize', 16)
legend('Senza Coda', 'Con Coda', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CL', 'FontSize', 12)
grid on
hold off

% Plot CDi-alpha
figure(4)
plot(alpha, CDi_C172, '-o', 'LineWidth', 1.2)
hold on
plot(alpha, CDi_C172_tail, '-o', 'LineWidth', 1.2)
title(sprintf('Curva CDi(%c) Cessna 172', char(0x03B1)), 'FontSize', 16)
legend('Senza Coda', 'Con Coda', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CDi', 'FontSize', 12)
grid on
hold off

% Plot Polare
%figure(5)
%plot(CDi_C172, CL_C172, '-o', 'LineWidth', 1.2)
%hold on
%plot(CDi_C172_tail, CL_C172_tail, '-o', 'LineWidth', 1.2)
%title('Polare Cessna 172', 'FontSize', 16)
%legend('Senza Coda', 'Con Coda', 'FontSize', 12)
%xlabel('CDi', 'FontSize', 12)
%ylabel('CL', 'FontSize', 12)
%grid on
%hold off

%% Confronto Tecnam P2008 (Con e Senza Piani di Coda)
% Plot CL-alpha
figure(3)
plot(alpha, CL_P2008, '-o', 'LineWidth', 2)
hold on
plot(alpha, CL_P2008_tail, '-o', 'LineWidth', 2)
title(sprintf('Curva CL(%c) Tecnam P2008', char(0x03B1)), 'FontSize', 16)
legend('Senza Coda', 'Con Coda', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CL', 'FontSize', 12)
grid on
hold off

% Plot CDi-alpha
figure(4)
plot(alpha, CDi_P2008, '-o', 'LineWidth', 2)
hold on
plot(alpha, CDi_P2008_tail, '-o', 'LineWidth', 2)
title(sprintf('Curva CDi(%c) Tecnam P2008', char(0x03B1)), 'FontSize', 16)
legend('Senza Coda', 'Con Coda', 'FontSize', 12)
xlabel(sprintf('%c (deg)', char(0x03B1)), 'FontSize', 12)
ylabel('CDi', 'FontSize', 12)
grid on
hold off

% Plot Polare
%figure(5)
%plot(CDi_P2008, CL_P2008, '-o', 'LineWidth', 1.2)
%hold on
%plot(CDi_P2008_tail, CL_P2008_tail, '-o', 'LineWidth', 1.2)
%title('Polare Tecnam 2008', 'FontSize', 16)
%legend('Senza Coda', 'Con Coda', 'FontSize', 12)
%xlabel('CDi', 'FontSize', 12)
%ylabel('CL', 'FontSize', 12)
%grid on
%hold off