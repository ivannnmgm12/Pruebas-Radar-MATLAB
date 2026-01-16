% 1. DEFINICION DE CTES

c = 3e8;  % vel luz

% 2. GENERACION PULSO RADAR

fc = 3e9;   % frecuencia portadora (Hz)
fs = 100e6; % frecuecnia de muestreo (Hz)
Tp = 0.2e-6    % ancho del pulso (s)

t = 0:1/fs:200e-6;  % ventana suficintemente grande para ver el eco

tx_pulse = rectpuls(t - Tp/2, Tp); % pulso rectangular

% Visualización pulso rectangular (pulso corto y limpio--> mejor resolucion
% en rango)

figure
plot(t*1e6, tx_pulse)      
xlabel('Tiempo (µs)')
ylabel('Amplitud')
title('Pulso radar transmitido') 
grid on

%3. SIMULACION DEL BLANCO

R_targets = [15e3 15.2e3 15.4e3];  % rango del blanco

rx_echo = zeros(size(t));  % inicializamos señal recibida

for k = 1:length(R_targets) % generar ecos individuales
    tau_k = 2 * R_targets(k) /c;
    rx_echo = rx_echo + rectpuls(t - tau_k - Tp/2, Tp);
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

SNR_dB = 5;  %Ruido termico. Radar con SNR bajo (realista)
rx_signal = awgn(rx_echo, SNR_dB, 'measured');

% visualizacion

figure
plot(t*1e6, rx_signal)
xlabel('Tiempo (µs)')
ylabel('Amplitud')
title('Señal recibida con múltiples blancos')
grid on

% 5. MATCHED FILTER--> Maximiza SNR

matched_filter = fliplr(tx_pulse);  % invierte el pulso en el tiempo
mf_output = conv(rx_signal, matched_filter, 'same');  % "busca" la forma del pulso en el ruido   

% visualizacion

figure
plot(t*1e6, abs(mf_output))
xlabel('Tiempo (µs)')
ylabel('|Salida MF|')
title('Matched filter - resolución en rango')
grid on

% 7. ESTIMACION DEL RANGO

[~, idx] = max(abs(mf_output));  % encontrar el pico
tau_est = t(idx); 

R_est = c * tau_est / 2;         % calcular rango estimado

fprintf('Rango estimado: %.2f km\n', R_est/1e3);  % resultado
