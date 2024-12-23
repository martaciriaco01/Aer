clc;
clear all;
close all;


%% carico i dati

alpha_range = -6:1:8;

load('file_cessna.mat');
load('file_tecnam.mat');



%% CL/alpha
figure;
plot(alpha_range, CL_values_Cessna,'-o', 'LineWidth', 1.2);
hold on;
plot(alpha_range, CL_values_Tecnam,'-o', 'LineWidth', 1.2);
xlabel('Angolo di Attacco (gradi)');
ylabel('C_L');
title('CL vs Angolo di Attacco', 'FontSize', 16);
grid on;
legend('Cessna 172 Skyhawk', 'Tecnam P2008', 'Location', 'Best');



%% Attrito indotto
figure;
plot(alpha_range, CD_values_Cessna, '-o', 'LineWidth', 1.2);
hold on;
plot(alpha_range, CD_values_Tecnam, '-o', 'LineWidth', 1.2);
xlabel('Angolo di Attacco (gradi)');
ylabel('C_{Di}');
title('CDi vs Angolo di Attacco', 'FontSize', 16);
grid on;
legend('Cessna 172 Skyhawk', 'Tecnam P2008', 'Location', 'Best');



%% Distribuzione circolazione rispetto all'ellittica
figure;
plot(y_Cessna, GammaDistribution_Cessna, 'o-',  'Color', [0 0.4470 0.7410], 'LineWidth', 1.2);
hold on;
plot(y_Cessna, Gamma_elliptic_Cessna, 'r--', 'LineWidth', 1.2);
hold on;
plot(y_Cessna, GammaDistribution_Tecnam, 'o-',  'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
hold on;
plot(y_Cessna, Gamma_elliptic_Tecnam, 'r--', 'LineWidth', 1.2);

title('Circulation Distribution over the Span', 'FontSize', 16);
xlabel('Spanwise Position [m]');
ylabel('Circulation (\Gamma) [m^2/s]');
legend('Computed Circulation', 'Elliptic Circulation');
grid on;
hold off;
legend('Cessna 172 Skyhawk', 'Circolazione ellittica', 'Tecnam P2008', 'Location', 'Best');




%% Polari
figure;
plot(CD_values_Cessna, CL_values_Cessna, '-o', 'LineWidth', 1.2);
hold on;
plot(CD_values_Tecnam, CL_values_Tecnam, '-o', 'LineWidth', 1.2);
xlabel('C_D');
ylabel('C_L');
title('Polare', 'FontSize', 16);
grid on;
legend('Cessna 172 Skyhawk', 'Tecnam P2008', 'Location', 'Best');
hold off;