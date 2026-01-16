% 1. DEFINICION DE CTES

c = 3e8;  % vel luz

% 2. GENERACION PULSO RADAR

Tp = 10e-6;        % Pulso LARGO (10 µs)
B  = 20e6;         % Ancho de banda (20 MHz)
fs = 100e6;        % Frecuencia de muestreo

t = 0:1/fs:200e-6;  % ventana temporal

tx_pulse = chirp(t, -B/2, Tp, B/2); % generacion del chirp transmitido
tx_pulse(t > Tp) = 0;   % ventana temporal

% Visualizacion TIEMPO-FRECUENCIA
spectrogram(tx_pulse, 256, 200, 256, fs, 'yaxis');
ylim([0 30])   % MHz
xlim([0 15])   % µs
title('Chirp transmitido (tiempo-frecuencia');

% Visualización chirp 

figure
plot(t*1e6, tx_pulse)      
xlabel('Tiempo (µs)')
ylabel('Amplitud')
title('Pulso radar transmitido') 
grid on

%3. SIMULACION DEL BLANCO

R_targets = [15e3 15.1e3];  % rango del blanco 

rx_echo = zeros(size(t));          % inicializamos señal recibida

for k = 1:length(R_targets)        % generamos multiples blancos
    tau_k = 2 * R_targets(k) / c;
    rx_echo = rx_echo + chirp(t - tau_k, -B/2, Tp, B/2) .* (t > tau_k & t < tau_k + Tp);
end

  
% comparacion visual

figure
plot(t*1e6, tx_pulse, 'b', 'LineWidth', 1.2)
hold on
plot(t*1e6, rx_echo, 'r', 'LineWidth', 1.2)
xlabel('Tiempo (µs)')
ylabel('Amplitud')
legend('Transmitido','Eco recibido')
title('Pulso transmitido y eco')
grid on

% 4. AÑADIMOS RUIDO

SNR_dB = 5;  % Ruido termico. Radar con SNR bajo (realista)
rx_signal = awgn(rx_echo, SNR_dB, 'measured');

% visualizacion

figure
plot(t*1e6, rx_signal)
xlabel('Tiempo (µs)')
ylabel('Amplitud')
title('Señal recibida con múltiples blancos')
grid on

% 5. MATCHED FILTER PARA CHIRP--> Maximiza SNR

matched_filter = fliplr(conj(tx_pulse));     % invierte el pulso en el tiempo
mf_output = conv(rx_signal, matched_filter, 'same');  % "busca" la forma del pulso en el ruido   

% visualizacion

figure
plot(t*1e6, abs(mf_output))
xlim([105 120])
xlabel('Tiempo (µs)')
ylabel('|Salida MF|')
title('Pulse compression con chirp')
grid on

% 7. ESTIMACION DEL RANGO

[~, idx] = max(abs(mf_output));  % encontrar el pico
tau_est = t(idx); 

R_est = c * tau_est / 2;         % calcular rango estimado

fprintf('Rango estimado: %.2f km\n', R_est/1e3);  % resultado
