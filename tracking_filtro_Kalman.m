%% ============================================================
%  EJERCICIO 5 - TRACKING DE BLANCO CON FILTRO DE KALMAN
%  Radar coherente: rango + velocidad
% ============================================================

clear
clc
close all

%% ============================================================
% 1. PARAMETROS DEL RADAR
% ============================================================

PRF = 1e3;                 % Frecuencia de repeticion (Hz)
T = 1/PRF;                 % Periodo entre pulsos
Np = 50;                   % Numero de medidas

%% ============================================================
% 2. BLANCO REAL (VERDAD TERRENO)
% ============================================================

r0 = 15e3;                 % Rango inicial (m)
v0 = 30;                   % Velocidad constante (m/s)

t = (0:Np-1) * T;

r_true = r0 + v0 * t;
v_true = v0 * ones(size(t));

%% ============================================================
% 3. MEDIDAS RADAR (CON RUIDO)
% ============================================================

sigma_r = 30;              % Desviacion error rango (m)
sigma_v = 1;               % Desviacion error velocidad (m/s)

r_meas = r_true + sigma_r * randn(size(t));
v_meas = v_true + sigma_v * randn(size(t));

z = [r_meas; v_meas];      % Vector de medidas

%% ============================================================
% 4. FILTRO DE KALMAN
% ============================================================

% Modelo dinamico
F = [1 T;
     0 1];

% Modelo de medida
H = eye(2);

% Covarianzas
Q = [1 0;
     0 0.5];               % Ruido de modelo

R = [sigma_r^2 0;
     0 sigma_v^2];         % Ruido de medida

% Inicializacion
x_est = [r_meas(1); v_meas(1)];
P = diag([100^2 10^2]);

x_history = zeros(2, Np);

%% ============================================================
% 5. BUCLE DE TRACKING
% ============================================================

for k = 1:Np

    % --- Prediccion ---
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;

    % --- Actualizacion ---
    K = P_pred * H' / (H * P_pred * H' + R);

    x_est = x_pred + K * (z(:,k) - H * x_pred);
    P = (eye(2) - K * H) * P_pred;

    x_history(:,k) = x_est;
end

%% ============================================================
% 6. VISUALIZACION - RANGO
% ============================================================

figure
plot(t, r_true/1e3, 'y--', 'LineWidth', 1.5)
hold on
plot(t, r_meas/1e3, 'rx')
plot(t, x_history(1,:)/1e3, 'b', 'LineWidth', 2)
xlabel('Tiempo (s)')
ylabel('Rango (km)')
legend('Real','Medido','Kalman')
title('Tracking de rango')
grid on

%% ============================================================
% 7. VISUALIZACION - VELOCIDAD
% ============================================================

figure
plot(t, v_true, 'y--', 'LineWidth', 1.5)
hold on
plot(t, v_meas, 'rx')
plot(t, x_history(2,:), 'b', 'LineWidth', 2)
xlabel('Tiempo (s)')
ylabel('Velocidad (m/s)')
legend('Real','Medido','Kalman')
title('Tracking de velocidad')
grid on

