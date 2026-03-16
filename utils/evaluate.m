function [LCOE_penalized, CF, AEP, cable] = evaluate(wf, turbine, cable, layout_coords)

pos = layout_coords';
T = turbine.turbine_num;

power_accum = zeros(T,1);

for it = 1:length(wf.theta)
    th = wf.theta(it);

    R = [cos(th) -sin(th); sin(th) cos(th)];
    rot = R * pos;

    for iv = 1:length(wf.velocity)

        v = wf.velocity(iv);
        ws = jensen_model(rot, turbine, v);
        p  = layout_power(ws, turbine);

        w = wf.f_theta_v(it,iv);
        power_accum = power_accum + p * w;

    end
end

AEP_t = power_accum * 365 * 24;
AEP = sum(AEP_t);
CF = AEP / (365*24*turbine.Pinst*T);

H = turbine.hub_height;
D = turbine.rotor_diameter;
SD = wf.sea_depth;

PWT = turbine.Pinst / 1000;
P_inst = PWT * T;

C_W  = 2.95e3 * log(PWT) - 375.2;
C_ist = 0.113 * C_W;
C_f = 320 * PWT * (1 + 0.02*(SD-8)) * (1 + 0.8e-6*(H*(D/2)^2 - 1e5));

toC = (C_W + C_ist + 1.5*C_f) * T;

if isfield(cable,'total_cost')
    C_C_inner = cable.total_cost;
else
    n_types = length(wf.innercable_price);
    len = cable.length_by_type;

    if numel(len) < n_types
        len = [len zeros(1,n_types-numel(len))];
    end

    C_C_inner = sum(wf.innercable_price .* len);
end

C_C_export = wf.offshore_length * wf.exportercable_price;
C_C = (C_C_inner + C_C_export) / 1000;

C_S_off = 539 * P_inst^0.678;
C_S_on  = 87.250 * P_inst;
C_S = C_S_off + C_S_on;

esC = C_S + C_C;

cost_per_mw = [7 92 88 52 144 130 89 109 365];
power = [495 120 108 60 NaN 40 40 NaN NaN];
years = [2004 2010 2007 2006 2010 2010 2005 2009 2010];

current_year = 2024;

yd = current_year - years;
valid = ~isnan(power);

w = power(valid) ./ yd(valid);

C_base = sum(cost_per_mw(valid).*w) / sum(w);
dpC = C_base * P_inst;

deC = 0.0093 * C_W * T;

C_O = 78.2;
C_O_y = C_O * P_inst;

Y = 25;

othC = 0.114*(toC+esC) + 0.092*C_O_y*Y;

r = 0.05;

num = 0;
den = 0;

for y = 1:Y
    num = num + C_O_y/(1+r)^y;
    den = den + AEP/(1+r)^y;
end

C_I = toC + esC + dpC + deC + othC;
num = num + C_I;

LCOE = num/den * 1e3;

LCOE_penalized = LCOE;

end


function power = layout_power(v, turbine)

power = zeros(size(v));

for i = 1:length(v)

    if v(i) < turbine.CutIn || v(i) > turbine.CutOut
        power(i) = 0;
    else
        power(i) = spline(turbine.ws, turbine.pc, v(i));
    end

end

end


function ws = jensen_model(pos, turbine, U0)

T = size(pos,2);

Ct = 0.8;
alpha = 0.1;
R = turbine.rotor_radius;

ws = U0 * ones(T,1);

[~,order] = sort(pos(2,:),'descend');

for j = 2:T

    tj = order(j);

    for i = 1:j-1

        ti = order(i);

        dx = pos(1,tj) - pos(1,ti);
        dy = pos(2,tj) - pos(2,ti);

        if dy < 0

            r = R + alpha*abs(dy);

            if abs(dx) < r
                decay = 1 - 2*a_func(Ct)*(R/r)^2;
                ws(tj) = ws(tj) * decay;
            end

        end

    end

end

end


function a = a_func(Ct)

a = 0.5 * (1 - sqrt(1 - Ct));

end