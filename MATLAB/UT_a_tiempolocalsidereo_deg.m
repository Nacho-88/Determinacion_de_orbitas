function lst = UT_a_tiempolocalsidereo_deg(dateTimeUTC, longitude)
    % Calcula el Tiempo Local Sidéreo (LST) dado una fecha y hora en UTC y la longitud geográfica.
    % Esta función devuelve el LST en GRADOS (deg).
    %
    % Parámetros:
    %   dateTimeUTC - Fecha y hora en UTC (tipo datetime).
    %   longitude   - Longitud geográfica del observador en grados 
    %                 (positivo hacia el este, negativo hacia el oeste).
    %
    % Retorno:
    %   lst - Tiempo Local Sidéreo (en horas).
    
    % Convertir la fecha y hora a formato Julian Date (JD)
    JD = juliandate(dateTimeUTC);
    
    % Calcular el número de días desde J2000.0
    T = (JD - 2451545.0) / 36525;  % Años julianos desde J2000.0
    
    % Calcular el tiempo sidéreo en Greenwich (GST) en grados
    GST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) ...
          + 0.000387933 * T^2 - T^3 / 38710000;
    
    % Ajustar GST al rango [0, 360) grados
    GST = mod(GST, 360);
    if GST < 0
        GST = GST + 360;
    end
    
    % Calcular el tiempo local sidéreo (LST) en grados
    LST = GST + longitude;
    
    % Ajustar LST al rango [0, 360) grados
    LST = mod(LST, 360);
    if LST < 0
        LST = LST + 360;
    end
    
    lst = LST;
end