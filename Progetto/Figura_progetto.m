figure();
hold on;

risultato_simulazione = sim("Simulink_progetto_non_linearizzato.slx");

angoli = risultato_simulazione.simout.Data(1:40:end);

% Rettangolo 1
x1 = [20 22 22 20];
y1 = [0 0 12 12];
patch(x1, y1, 'b', 'FaceAlpha', 1);


% Rettangolo 2
x2 = [33 35 35 33];
y2 = [0 0 12 12];
patch(x2, y2, 'b', 'FaceAlpha', 1);

% Rettangolo 3
x3 = [19 36 36 19];
y3 = [10 10 12 12];
patch(x3, y3, 'b', 'FaceAlpha', 1);

% Rettangolo 4
x4 = [8 19 19 8];
y4 = [15 10 12 17];
patch(x4, y4, 'b', 'FaceAlpha', 1);

% Parametri della circonferenza
center_x = 9;
center_y = 13;
radius = 5;

% Creazione della circonferenza
theta_circle = linspace(0, 2*pi, 100);
x_circle = center_x + radius * cos(theta_circle);
y_circle = center_y + radius * sin(theta_circle);
plot(x_circle, y_circle, 'k', 'LineWidth', 1.5);

% Parametri della spazzata
sweep_color = [1 0 0];  % Colore dell'area spazzata

% Punto che ruota sulla circonferenza
point_radius = 5;
theta_point = 0;  % Fase iniziale

% Inizializza la matrice la posizione del punto
sweep_trace = [];

% Creazione del punto e della traccia dell'area spazzata
point_x = center_x + point_radius * cos(theta_point);
point_y = center_y + point_radius * sin(theta_point);
point_handle = plot(point_x, point_y, 'ro', 'MarkerSize', 10);

% Impostazioni dell'asse
axis equal;
axis([0 40 0 20]);

xlabel('X');
ylabel('Y');

title('Disegno con Patch in MATLAB');

% Legenda
legend('Rettangolo 1', 'Rettangolo 2', 'Rettangolo 3', 'Rettangolo 4', 'Circonferenza', 'Punto', 'Location', 'Best', 'AutoUpdate', 'off');
grid on;

% Animazione della rotazione del punto sulla circonferenza
for i = 1:length(angoli)
    theta_point = pi/2 - angoli(i); % Imposta l'angolo in base ai dati di simulazione
    point_x = center_x + point_radius * cos(theta_point);
    point_y = center_y + point_radius * sin(theta_point);
    set(point_handle, 'XData', point_x, 'YData', point_y);
    
    % Aggiungi il segmento dell'area spazzata sulla circonferenza
    sweep_trace = [sweep_trace; point_x, point_y];
    plot(sweep_trace(:, 1), sweep_trace(:, 2), 'Color', sweep_color, 'LineWidth', 1.5);
    
    pause(0.1);
end
hold off;