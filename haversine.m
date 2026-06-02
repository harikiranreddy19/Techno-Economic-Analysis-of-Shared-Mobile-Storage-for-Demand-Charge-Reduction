function distance = haversine(lat1, lon1, lat2, lon2)
    % obtaining distance in miles

    R = 6371.0; % Radius of the Earth in kilometers
    
    % Convert latitude and longitude from degrees to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);
    
    % Differences in coordinates
    dlon = lon2 - lon1;
    dlat = lat2 - lat1;
    
    % Haversine formula
    a = sin(dlat / 2).^2 + cos(lat1) .* cos(lat2) .* sin(dlon / 2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    % Calculate the distance in miles
    sign = 1*double(lat1<lat2) - 1*double(lat1>lat2) + ...
           1*double(lon1<lon2) - 1*double(lon1>lon2);

    distance = sign.*0.6214 *R .* c;
end