function radians = hms_to_rad(hours, minutes, seconds)

    % Convierte un ángulo en horas, minutos y segundos a radianes
    % 
    % Parámetros:
    %   hours   - Número de horas (entero o decimal)
    %   minutes - Número de minutos (entero o decimal)
    %   seconds - Número de segundos (entero o decimal)
    %
    % Retorno:
    %   radians - Ángulo convertido a radianes
    
    % Conversión a grados
    degrees = (hours * 15) + (minutes * 15 / 60) + (seconds * 15 / 3600);
    
    % Conversión a radianes
    radians = deg2rad(degrees);
end