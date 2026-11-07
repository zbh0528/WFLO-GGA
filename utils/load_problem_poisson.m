function [wf, turbine, wind_name, turbine_name] = load_problem_poisson(case_name)

switch case_name
    case 'China_Shanghai_Lingang'
        wind_file = './data/wind/China_Shanghai_Lingang.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'China_Zhuhai_Guishan_Hai'
        wind_file = './data/wind/China_Zhuhai_Guishan_Hai.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Denmark_Horns_Rev_1'
        wind_file = './data/wind/Denmark_Horns_Rev_1.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Denmark_Nysted'
        wind_file = './data/wind/Denmark_Nysted.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Netherlands_Egmond_aan_Zee'
        wind_file = './data/wind/Netherlands_Egmond_aan_Zee.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Netherlands_Prinses_Amaliawindpark'
        wind_file = './data/wind/Netherlands_Prinses_Amaliawindpark.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'UK_London_Array'
        wind_file = './data/wind/UK_London_Array.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'UK_Sheringham_Shoal'
        wind_file = './data/wind/UK_Sheringham_Shoal.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Denmark_Rodsand_II'
        wind_file = './data/wind/Denmark_Rodsand_II.mat';
        sea_depth = 7; offshore_length = 0; alpha = 0.3;

    case 'Denmark-Anholt'
        wind_file = './data/wind/windlib/gwc_denmark_anholt.mat';
        sea_depth = 32; offshore_length = 6000; alpha = 1;

    otherwise
        error('Unknown case');
end


turbine_file = './data/turbine/Vetas_4200kwturbine_speed-pv-Ct_value.csv';

[~, wind_name] = fileparts(wind_file);
[~, turbine_name] = fileparts(turbine_file);

T = readtable(turbine_file,'Range','A2:C46');

turbine.ws = table2array(T(:,1));
turbine.pc = table2array(T(:,2));
turbine.Ct = table2array(T(:,3));

turbine.Pinst = 4200;
turbine.hub_height = 91.5;
turbine.rotor_diameter = 117.0;
turbine.rotor_radius = turbine.rotor_diameter/2;
turbine.CutIn = 3;
turbine.CutOut = 25;
turbine.turbulence_intensity = 0.1;


layout0 = load_layout([case_name '.csv']) * 1000;
turbine.turbine_num = size(layout0,1);

wf.layout0 = layout0;

wind = load(wind_file);
wf.velocity = wind.velocity;
wf.f_theta_v = wind.f_theta_v';
wf.theta = wind.theta;

wf.name = case_name;
wf.surface_roughness = 0.25e-3;
wf.sea_depth = sea_depth;
wf.offshore_length = offshore_length;


K = boundary(layout0(:,1),layout0(:,2),alpha);

bx = layout0(K,1);
by = layout0(K,2);

boundary_polygon = [bx by];


r = 3 * turbine.rotor_diameter;

bbox = [min(bx) max(bx); min(by) max(by)];

while true

    pts = bridson_poisson_sampler(r,bbox,boundary_polygon,30);

    if size(pts,1) >= turbine.turbine_num
        break
    end

    r = r * 0.99;

end


wf.min_turbine_dis = r;
wf.candidate_points = pts;
wf.N_candidate = size(pts,1);

wf.boundary = boundary_polygon;
wf.substation = mean(layout0,1);


wf.innercable_voltage_kV = 33;
wf.exportercable_voltage_kV = 220;
wf.cos_phi = 0.95;

wf.turbine_power_MW = turbine.Pinst/1000;

I_turbine = wf.turbine_power_MW*1e6/(sqrt(3)*wf.innercable_voltage_kV*1e3*wf.cos_phi);

wf.innercable_section_area = [150 400 800];
wf.innercable_rated_current = [430 680 900];
wf.innercable_price = [192 321 506];

wf.cable_capacity = floor(wf.innercable_rated_current / I_turbine);

wf.exportercable_section_area = 1600;
wf.exportercable_price = 601.500;

end



function pts = bridson_poisson_sampler(r,bbox,boundary,k)

w = r/sqrt(2);

cols = ceil((bbox(1,2)-bbox(1,1))/w);
rows = ceil((bbox(2,2)-bbox(2,1))/w);

grid = -ones(rows,cols);

pts = [];
active = [];

while true

    p0 = [rand*(bbox(1,2)-bbox(1,1))+bbox(1,1) ...
          rand*(bbox(2,2)-bbox(2,1))+bbox(2,1)];

    if inpolygon(p0(1),p0(2),boundary(:,1),boundary(:,2))

        pts = p0;
        active = 1;

        row = floor((p0(2)-bbox(2,1))/w)+1;
        col = floor((p0(1)-bbox(1,1))/w)+1;

        grid(row,col) = 1;

        break

    end
end


while ~isempty(active)

    idx = active(randi(length(active)));
    base = pts(idx,:);

    found = false;

    for i = 1:k

        theta = 2*pi*rand;
        m = r*(1+rand);

        pt = base + m*[cos(theta) sin(theta)];

        if pt(1)<bbox(1,1) || pt(1)>bbox(1,2) || pt(2)<bbox(2,1) || pt(2)>bbox(2,2)
            continue
        end

        if ~inpolygon(pt(1),pt(2),boundary(:,1),boundary(:,2))
            continue
        end

        row = floor((pt(2)-bbox(2,1))/w)+1;
        col = floor((pt(1)-bbox(1,1))/w)+1;

        if row<1 || col<1 || row>rows || col>cols
            continue
        end

        neigh = [];

        for ii=-2:2
            for jj=-2:2

                rr = row+ii;
                cc = col+jj;

                if rr>=1 && cc>=1 && rr<=rows && cc<=cols

                    id = grid(rr,cc);

                    if id>0
                        neigh(end+1,:) = pts(id,:);
                    end

                end
            end
        end

        if isempty(neigh) || all(vecnorm(neigh-pt,2,2)>=r)

            pts(end+1,:) = pt;
            active(end+1) = size(pts,1);

            grid(row,col) = size(pts,1);

            found = true;
            break

        end
    end

    if ~found
        active(active==idx) = [];
    end

end

end