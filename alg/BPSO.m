function best_fit = BPSO(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <1997> <single> <binary> <none>
% Binary particle swarm optimization
% w  --- 0.7 --- Inertia weight controlling velocity memory
% c1 --- 1.5 --- Cognitive acceleration coefficient
% c2 --- 1.5 --- Social acceleration coefficient
%
%------------------------------- Reference --------------------------------
% J. Kennedy and R. C. Eberhart. A discrete binary version of the particle
% swarm algorithm. Proceedings of the IEEE International Conference on
% Systems, Man, and Cybernetics, 1997: 4104-4108.
%
% BPSO: binary particle swarm optimization with fixed Top-M site selection
% - Chromosome representation: a length-N binary vector with exactly M ones
% - Velocity update: standard PSO-style update
% - Position update: choose the top-M velocity entries
% - Supports unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

%% ======================== Parameter preparation ========================
M = turbine.turbine_num;
N = wf.N_candidate;

% PSO parameters
w  = 0.7;
c1 = 1.5;
c2 = 1.5;

if nargin < 9 || isempty(routing_fn)
    error('BPSO:MissingRoutingFunction', ...
        'routing_fn must be provided, e.g., @cr_mst, @cr_sweep, or @cr_sector.');
end

%% ======================== Output folder ========================
mat_dir = results_dir;
if ~exist(mat_dir, 'dir')
    mkdir(mat_dir);
end

%% ======================== 1) Population initialization ========================
particles = repmat(struct('Position', [], 'Velocity', []), popsize, 1);
pBest     = repmat(struct('Position', [], 'Cost', inf), popsize, 1);
gBest     = struct('Position', [], 'Cost', inf, 'AEP', [], 'CF', [], 'CableInfo', []);

% 1.1 Initial positions
if nargin >= 8 && ~isempty(Pop0)
    if size(Pop0, 2) == M
        for i = 1:popsize
            pos = zeros(1, N);
            idx = Pop0(i, :);
            pos(idx) = 1;
            pos = ensure_cardinality_M(pos, M);
            particles(i).Position = pos;
        end
    elseif size(Pop0, 2) == N
        for i = 1:popsize
            pos = Pop0(i, :);
            pos = ensure_cardinality_M(pos, M);
            particles(i).Position = pos;
        end
    else
        error('Pop0 size mismatch: expected popsize-by-M or popsize-by-N.');
    end
else
    for i = 1:popsize
        pos = zeros(1, N);
        idx = randperm(N, M);
        pos(idx) = 1;
        particles(i).Position = pos;
    end
end

% 1.2 Initial velocities
for i = 1:popsize
    particles(i).Velocity = zeros(1, N);
end

%% ======================== 2) Initial evaluation ========================
fprintf('Initializing BPSO population evaluation...\n');

fit_init    = zeros(popsize, 1);
CF_init     = zeros(popsize, 1);
AEP_init    = zeros(popsize, 1);
cables_init = cell(popsize, 1);
routing_time_init = 0;

for i = 1:popsize
    ind = find(particles(i).Position == 1);
    [cost, cf, aep, cable, rt] = evaluate_individual(ind, wf, turbine, routing_fn);
    routing_time_init = routing_time_init + rt;

    pBest(i).Position = particles(i).Position;
    pBest(i).Cost     = cost;

    if cost < gBest.Cost
        gBest.Cost      = cost;
        gBest.Position  = particles(i).Position;
        gBest.AEP       = aep;
        gBest.CF        = cf;
        gBest.CableInfo = cable;
    end

    fit_init(i)    = cost;
    CF_init(i)     = cf;
    AEP_init(i)    = aep;
    cables_init{i} = cable;

    fprintf('Initial individual %2d: fitness = %.6f\n', i, cost);
end

%% Save generation 1
generations = cell(max_it, 1);
population_idx = positions_to_indices_matrix(particles, M);

gen_struct = struct();
gen_struct.population        = population_idx;
gen_struct.fitness           = fit_init;
gen_struct.CF                = CF_init;
gen_struct.AEP               = AEP_init;
gen_struct.cables            = cables_init;
gen_struct.runtime           = 0;
gen_struct.runtime_routing   = routing_time_init;
gen_struct.runtime_evolution = 0;
generations{1} = gen_struct;

best_fit = gBest.Cost;
fprintf('[BPSO] Generation %3d: best = %.8f | routing time = %.2fs\n', best_fit ~= inf, best_fit, routing_time_init);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting BPSO evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    fit_gen    = zeros(popsize, 1);
    CF_gen     = zeros(popsize, 1);
    AEP_gen    = zeros(popsize, 1);
    cables_gen = cell(popsize, 1);

    for i = 1:popsize
        r1 = rand(1, N);
        r2 = rand(1, N);

        particles(i).Velocity = ...
            w  * particles(i).Velocity + ...
            c1 * r1 .* (pBest(i).Position - particles(i).Position) + ...
            c2 * r2 .* (gBest.Position   - particles(i).Position);

        [~, sort_idx] = sort(particles(i).Velocity, 'descend');
        new_pos = zeros(1, N);
        new_pos(sort_idx(1:M)) = 1;

        ind = find(new_pos == 1);
        [cost, cf, aep, cable, rt] = evaluate_individual(ind, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;

        fit_gen(i)    = cost;
        CF_gen(i)     = cf;
        AEP_gen(i)    = aep;
        cables_gen{i} = cable;

        if cost < pBest(i).Cost
            pBest(i).Cost     = cost;
            pBest(i).Position = new_pos;
        end

        if cost < gBest.Cost
            gBest.Cost      = cost;
            gBest.Position  = new_pos;
            gBest.AEP       = aep;
            gBest.CF        = cf;
            gBest.CableInfo = cable;
        end

        particles(i).Position = new_pos;
    end

    total_time     = toc(gen_tic);
    evolution_time = total_time - routing_time_this_gen;

    population_idx = positions_to_indices_matrix(particles, M);

    gen_struct = struct();
    gen_struct.population        = population_idx;
    gen_struct.fitness           = fit_gen;
    gen_struct.CF                = CF_gen;
    gen_struct.AEP               = AEP_gen;
    gen_struct.cables            = cables_gen;
    gen_struct.runtime           = total_time;
    gen_struct.runtime_routing   = routing_time_this_gen;
    gen_struct.runtime_evolution = evolution_time;
    generations{gen} = gen_struct;

    best_fit = gBest.Cost;
    fprintf('[BPSO] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('BPSO finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function pos = ensure_cardinality_M(pos, M)
N = numel(pos);
k = sum(pos);
if k == M
    return;
end

idx1 = find(pos == 1);
idx0 = find(pos == 0);

if k > M
    turn_off = randsample(idx1, k - M, false);
    pos(turn_off) = 0;
else
    turn_on = randsample(idx0, M - k, false);
    pos(turn_on) = 1;
end

if sum(pos) ~= M
    pos = zeros(1, N);
    idx = randperm(N, M);
    pos(idx) = 1;
end
end

function pop_idx = positions_to_indices_matrix(particles, M)
popsize = numel(particles);
pop_idx = zeros(popsize, M);

for i = 1:popsize
    idx = find(particles(i).Position == 1);

    if numel(idx) > M
        idx = idx(1:M);
    elseif numel(idx) < M
        idx(end+1:M) = 1;
    end

    pop_idx(i, :) = sort(idx);
end
end

function [fit, cf, aep, cab, routing_time] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
layout_coords = wf.candidate_points(ind, :);
cable = routing_fn(layout_coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, layout_coords);
routing_time = toc(tic_routing);
end