function best_fit = AGPSO(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2022> <single> <integer> <none>
% Adaptive replacement strategy-incorporated particle swarm optimization
% c  --- 1.5  --- Acceleration coefficient controlling velocity update
% sg --- 7    --- Stagnation threshold triggering tournament-based revival
% pm --- 0.01 --- Mutation probability in exemplar-based crossover
%
%------------------------------- Reference --------------------------------
% Z. Lei, S. Gao, Y. Wang, Y. Yu, and L. Guo. An adaptive replacement
% strategy-incorporated particle swarm optimizer for wind farm layout
% optimization. Energy Conversion and Management, 2022, 269: 116174.
%
% AGPSO: adaptive genetic particle swarm optimization
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

M = turbine.turbine_num;
N = wf.N_candidate;

c = 1.5;

if nargin < 9 || isempty(routing_fn)
    error('AGPSO:MissingRoutingFunction', ...
        'routing_fn must be provided, e.g., @cr_mst, @cr_sweep, or @cr_sector.');
end

%% ======================== Output folder ========================
mat_dir = results_dir;
if ~exist(mat_dir, 'dir')
    mkdir(mat_dir);
end

%% ======================== Particle structures ========================
particles = repmat(struct('Position', [], 'Velocity', [], 'stagnate', 0), popsize, 1);
pBest     = repmat(struct('Position', [], 'Cost', inf), popsize, 1);
gBest     = struct('Position', [], 'Cost', inf, 'AEP', [], 'CF', [], 'CableInfo', []);

generations = cell(max_it, 1);

%% ======================== 1) Initialization ========================
fprintf('Initializing AGPSO population evaluation...\n');
routing_time_this_gen = 0;

fit_gen    = zeros(popsize, 1);
CF_gen     = zeros(popsize, 1);
AEP_gen    = zeros(popsize, 1);
cables_gen = cell(popsize, 1);

for i = 1:popsize
    if nargin >= 8 && ~isempty(Pop0)
        idx = Pop0(i, :);
    else
        idx = randperm(N, M);
    end

    particles(i).Position = unique_fix(idx, N);
    particles(i).Velocity = zeros(1, M);
    particles(i).stagnate = 0;

    [cost, cf, aep, cable, rt] = evaluate_individual(particles(i).Position, wf, turbine, routing_fn);
    routing_time_this_gen = routing_time_this_gen + rt;

    pBest(i).Position = particles(i).Position;
    pBest(i).Cost     = cost;

    if cost < gBest.Cost
        gBest.Cost      = cost;
        gBest.Position  = particles(i).Position;
        gBest.AEP       = aep;
        gBest.CF        = cf;
        gBest.CableInfo = cable;
    end

    fit_gen(i)    = cost;
    CF_gen(i)     = cf;
    AEP_gen(i)    = aep;
    cables_gen{i} = cable;

    fprintf('Initial individual %2d: fitness = %.6f\n', i, cost);
end

generation_struct = struct();
generation_struct.population        = vertcat(particles.Position);
generation_struct.fitness           = fit_gen;
generation_struct.CF                = CF_gen;
generation_struct.AEP               = AEP_gen;
generation_struct.cables            = cables_gen;
generation_struct.runtime           = 0;
generation_struct.runtime_routing   = routing_time_this_gen;
generation_struct.runtime_evolution = 0;
generations{1} = generation_struct;

best_fit = gBest.Cost;
fprintf('[AGPSO] Generation %3d: best = %.8f | routing = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

fprintf('\nStarting AGPSO evolution...\n');

%% ======================== AGPSO parameters ========================
sg = 7;
pm = 0.01;

%% ======================== 2) Evolution loop ========================
for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    fit_gen    = zeros(popsize, 1);
    CF_gen     = zeros(popsize, 1);
    AEP_gen    = zeros(popsize, 1);
    cables_gen = cell(popsize, 1);

    for i = 1:popsize
        % --- Exemplar-based crossover ---
        offsPbest = zeros(1, M);
        for d = 1:M
            k = randi([1, popsize]);
            if pBest(i).Cost < pBest(k).Cost
                r = rand;
                offsPbest(d) = round(r * pBest(i).Position(d) + (1 - r) * gBest.Position(d));
            else
                offsPbest(d) = pBest(k).Position(d);
            end

            if rand < pm
                offsPbest(d) = randi([1, N]);
            end
        end
        offsPbest = unique_fix(offsPbest, N);

        [cost_off, ~, ~, ~, rt_off] = evaluate_individual(offsPbest, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt_off;

        if cost_off < pBest(i).Cost
            pBest(i).Cost     = cost_off;
            pBest(i).Position = offsPbest;
            particles(i).stagnate = 0;
        else
            particles(i).stagnate = particles(i).stagnate + 1;
        end

        % --- Tournament-based revival ---
        if particles(i).stagnate >= sg
            tournament_size = max(2, round(0.2 * popsize));
            competitors = randperm(popsize, tournament_size);
            [~, winId] = min([pBest(competitors).Cost]);
            winner_idx = competitors(winId);

            pBest(i).Position = pBest(winner_idx).Position;
            pBest(i).Cost     = pBest(winner_idx).Cost;
            particles(i).stagnate = 0;
        end

        % --- Velocity update ---
        for d = 1:M
            particles(i).Velocity(d) = c * rand * (pBest(i).Position(d) - particles(i).Position(d));
        end

        % --- Position update ---
        new_pos = particles(i).Position + round(particles(i).Velocity);
        new_pos = unique_fix(new_pos, N);

        [cost, cf, aep, cab, rt] = evaluate_individual(new_pos, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;

        fit_gen(i)    = cost;
        CF_gen(i)     = cf;
        AEP_gen(i)    = aep;
        cables_gen{i} = cab;

        if cost < gBest.Cost
            gBest.Cost      = cost;
            gBest.Position  = new_pos;
            gBest.AEP       = aep;
            gBest.CF        = cf;
            gBest.CableInfo = cab;
        end

        particles(i).Position = new_pos;
    end

    total_time     = toc(gen_tic);
    evolution_time = total_time - routing_time_this_gen;

    generation_struct = struct();
    generation_struct.population        = vertcat(particles.Position);
    generation_struct.fitness           = fit_gen;
    generation_struct.CF                = CF_gen;
    generation_struct.AEP               = AEP_gen;
    generation_struct.cables            = cables_gen;
    generation_struct.runtime           = total_time;
    generation_struct.runtime_routing   = routing_time_this_gen;
    generation_struct.runtime_evolution = evolution_time;
    generations{gen} = generation_struct;

    fprintf('[AGPSO] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, gBest.Cost, total_time, evolution_time, routing_time_this_gen);

    best_fit = gBest.Cost;
end

%% ======================== 3) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('AGPSO finished. Best LCOE = %.10f\n', best_fit);
end

%% ======================== Helper functions ========================
function fixed = unique_fix(ind, N)
ind = max(1, min(N, round(ind)));
u = unique(ind, 'stable');
M = numel(ind);

if numel(u) == M
    fixed = u;
    return;
end

need = M - numel(u);
pool = setdiff(1:N, u);

if numel(pool) < need
    extra = setdiff(1:N, u);
    pool = unique([pool, extra], 'stable');
end

add = datasample(pool, need, 'Replace', false);
fixed = [u, add];
end

function [fit, cf, aep, cab, routing_time] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
layout_coords = wf.candidate_points(ind, :);
cable = routing_fn(layout_coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, layout_coords);
routing_time = toc(tic_routing);
end