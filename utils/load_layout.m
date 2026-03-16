function layout = load_layout(file)
data = readtable(file);

lat = data.centr_lat;
lon = data.centr_lon;

R = 6371; % km
lat0 = lat(1);
lon0 = lon(1);

x = (lon - lon0) * pi/180 * R * cosd(lat0);
y = (lat - lat0) * pi/180 * R;

layout = [x y];
end