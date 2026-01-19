%% ============================================================
%  RADAR CHIRP CON ESTIMACION CORRECTA DE DOPPLER
%  - Rango por matched filter
%  - Doppler por FFT en slow-time (celda detectada)
% ============================================================

clear
clc
close all

%% ============================================================
% 1. CONSTANTES FISICAS
% ============================================================

c = 3e8;                      % Velocidad de la luz (m/s)

%% ============================================================
% 2. PARAMETROS DEL RADAR
% ============================================================

fc = 3e9;                     % Frecuencia portadora (Hz)
lambda = c / fc;              % Longitud de onda (m)

Tp = 10e-6;                   % Duracion del chirp (10 us)
B  = 20e6;                    % Ancho de banda (20 MHz)
fs = 100e6;                   % Frecuencia de muestreo (Hz)

PRF = 1e3;                    % Frecuencia de repeticion de pulso (Hz)
Np  = 64;                     % Numero de pulsos coherentes

%% ============================================================
% 3. EJE TEMPORAL (FAST TIME)
% ============================================================

t = 0:1/fs:200e-6;            % Ventana temporal (~30 km)

%% ============================================================
% 4. PULSO CHIRP TRANSMITIDO
% ============================================================

tx_pulse = chirp(t, -B/2, Tp, B/2);
tx_pulse(t > Tp) = 0;

figure
spectrogram(tx_pulse, 256, 200, 256, fs, 'yaxis');
ylim([0 30])
xlim([0 15])
title('Chirp transmitido (tiempo–frecuencia)')

%% ============================================================
% 5. ESCENARIO DEL BLANCO
% ============================================================

R0 = 15e3;                    % Rango inicial (m)
v  = 30;                      % Velocidad radial (m/s)

%% ============================================================
% 6. GENERACION DE ECOS (Np x fast-time)
% ============================================================

rx_matrix = zeros(Np, length(t));

for p = 1:Np

    tp = (p-1) / PRF;                 % Tiempo lento
    Rp = R0 + v * tp;                 % Rango instantaneo
    tau_p = 2 * Rp / c;               % Retardo

    % Fase Doppler
    doppler_phase = exp(1j * 2*pi * (2*v/lambda) * tp);

    % Eco chirp retardado
    echo = chirp(t - tau_p, -B/2, Tp, B/2) .* ...
           (t > tau_p & t < tau_p + Tp);

    rx_matrix(p, :) = doppler_phase * echo;
end

% Ruido termico
rx_matrix = awgn(rx_matrix, 5, 'measured');

%% ============================================================
% 7. MATCHED FILTER (PROCESADO EN RANGO)
% ============================================================

mf = fliplr(conj(tx_pulse));

range_profiles = zeros(size(rx_matrix));

for p = 1:Np
    range_profiles(p,:) = conv(rx_matrix(p,:), mf, 'same');
end

%% ============================================================
% 8. EJE DE RANGO
% ============================================================

range_axis = c * t / 2;       % Rango fisico (m)

%% ============================================================
% 9. DETECCION DE LA CELDA DE RANGO DEL BLANCO
% ============================================================

% Integracion no coherente en slow-time
range_energy = mean(abs(range_profiles), 1);

[~, idx_r] = max(range_energy);

fprintf('Celda de rango detectada: %.2f km\n', range_axis(idx_r)/1e3);

%% ============================================================
% 10. EXTRACCION DE LA SEÑAL EN SLOW-TIME
% ============================================================

slow_time_signal = range_profiles(:, idx_r);

%% ============================================================
% 11. FFT DOPPLER (ESTIMACION CORRECTA)
% ============================================================

doppler_spectrum = fftshift(fft(slow_time_signal));

doppler_dB = 20*log10(abs(doppler_spectrum));
doppler_dB = doppler_dB - max(doppler_dB);   % Normalizacion

%% ============================================================
% 12. EJE DE VELOCIDAD
% ============================================================

fd = linspace(-PRF/2, PRF/2, Np);
v_axis = fd * lambda / 2;

%% ============================================================
% 13. VISUALIZACION ESPECTRO DOPPLER
% ============================================================

figure
plot(v_axis, doppler_dB, 'LineWidth', 1.5)
xlabel('Velocidad (m/s)')
ylabel('Nivel (dB)')
title('Espectro Doppler del blanco')
grid on
ylim([-40 0])

%% ============================================================
% 14. ESTIMACION DE VELOCIDAD
% ============================================================

[~, idx_v] = max(doppler_dB);

fprintf('Velocidad estimada: %.2f m/s\n', v_axis(idx_v));
