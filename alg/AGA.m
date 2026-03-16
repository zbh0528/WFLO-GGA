function best_fit = AGA(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2019> <single> <integer> <none>
% Adaptive genetic algorithm
% Cp --- 0.60 --- Crossover probability controlling recombination intensity
% Mp --- 0.10 --- Mutation probability applied to offspring individuals
% Ep --- 0.20 --- Proportion of elite individuals preserved each generation
% r  --- 0.50 --- Probability for selecting non-elite individuals into parent pool
%
%------------------------------- Reference --------------------------------
% X. Ju and F. Liu. Wind farm layout optimization using self-informed
% genetic algorithm with information guided exploitation. Applied Energy,
% 2019, 248: 429–445.
%
% AGA: adaptive genetic algorithm
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

%% ======================== Parameter preparation ========================
M = turbine.turbine_num;
N = wf.N_candidate;

% AGA parameters
Cp = 0.60; %#ok<NASGU>
Mp = 0.10;
Ep = 0.20;
r  = 0.50;
elite_num = max(1, floor(Ep * popsize));

if nargin < 9 || isempty(routing_fn)
    error('AGA:MissingRoutingFunction', ...
        'routing_fn must be provided, e.g., @cr_mst, @cr_sweep, or @cr_sector.');
end

%% ======================== Output folder ========================
mat_dir = results_dir;
if ~exist(mat_dir, 'dir')
    mkdir(mat_dir);
end

%% ======================== 1) Population initialization ========================
if nargin >= 8 && ~isempty(Pop0)
    pop = Pop0;
else
    pop = zeros(popsize, M);
    for i = 1:popsize
        pop(i, :) = randperm(N, M);
    end
end

%% ======================== 2) Initial evaluation ========================
fprintf('Initializing AGA population evaluation...\n');
fitness   = zeros(popsize, 1);
CF_array  = zeros(popsize, 1);
AEP_array = zeros(popsize, 1);
cables    = cell(popsize, 1);

routing_time_this_gen = 0;
for i = 1:popsize
    [fitness(i), CF_array(i), AEP_array(i), cables{i}, rt] = ...
        evaluate_individual(pop(i, :), wf, turbine, routing_fn);
    fprintf('Initial individual %2d: fitness = %.6f\n', i, fitness(i));
    routing_time_this_gen = routing_time_this_gen + rt;
end

%% Save generation 1
generations = cell(max_it, 1);
generation_struct = struct();
generation_struct.population        = pop;
generation_struct.fitness           = fitness;
generation_struct.CF                = CF_array;
generation_struct.AEP               = AEP_array;
generation_struct.cables            = cables;
generation_struct.runtime           = 0;
generation_struct.runtime_routing   = routing_time_this_gen;
generation_struct.runtime_evolution = 0;
generations{1} = generation_struct;

[best_fit, ~] = min(fitness);
fprintf('[AGA] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting AGA evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    % --- 3.1 Sort and extract elites ---
    [fitness, order] = sort(fitness, 'ascend');
    pop       = pop(order, :);
    CF_array  = CF_array(order);
    AEP_array = AEP_array(order);
    cables    = cables(order);

    elites_pop  = pop(1:elite_num, :);
    elites_fit  = fitness(1:elite_num);
    elites_CF   = CF_array(1:elite_num);
    elites_AEP  = AEP_array(1:elite_num);
    elites_cabs = cables(1:elite_num);

    % --- 3.2 Build parent pool ---
    parent_pool_idx = build_parent_pool_Aga(fitness, elite_num, r);
    parent_pool = pop(parent_pool_idx, :);
    if size(parent_pool, 1) < 2
        parent_pool = pop(1:min(2, popsize), :);
    end

    % --- 3.3 Generate next population ---
    new_pop    = zeros(popsize, M);
    new_fit    = zeros(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);

    % Directly copy elites
    new_pop(1:elite_num, :) = elites_pop;
    new_fit(1:elite_num)    = elites_fit;
    new_CF(1:elite_num)     = elites_CF;
    new_AEP(1:elite_num)    = elites_AEP;
    for e = 1:elite_num
        new_cables{e} = elites_cabs{e};
    end

    % Crossover + mutation
    idx_new = elite_num + 1;
    while idx_new <= popsize
        p1 = parent_pool(randi(size(parent_pool, 1)), :);
        p2 = parent_pool(randi(size(parent_pool, 1)), :);

        while isequal(p1, p2) && size(parent_pool, 1) > 1
            p2 = parent_pool(randi(size(parent_pool, 1)), :);
        end

        child = aga_single_point_crossover(p1, p2, N);

        if rand < Mp
            child = aga_point_mutation(child, N);
        end

        child = unique_fix(child, N);

        [new_fit(idx_new), new_CF(idx_new), new_AEP(idx_new), ...
            new_cables{idx_new}, rt] = evaluate_individual(child, wf, turbine, routing_fn);

        routing_time_this_gen = routing_time_this_gen + rt;
        new_pop(idx_new, :) = child;
        idx_new = idx_new + 1;
    end

    % Replace old population
    pop       = new_pop;
    fitness   = new_fit;
    CF_array  = new_CF;
    AEP_array = new_AEP;
    cables    = new_cables;

    % --- 3.4 Save this generation ---
    total_time     = toc(gen_tic);
    evolution_time = total_time - routing_time_this_gen;

    generation_struct.population        = pop;
    generation_struct.fitness           = fitness;
    generation_struct.CF                = CF_array;
    generation_struct.AEP               = AEP_array;
    generation_struct.cables            = cables;
    generation_struct.runtime           = total_time;
    generation_struct.runtime_routing   = routing_time_this_gen;
    generation_struct.runtime_evolution = evolution_time;
    generations{gen} = generation_struct;

    % --- 3.5 Logging ---
    [best_fit, ~] = min(fitness);
    fprintf('[AGA] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');
fprintf('AGA finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function parent_pool_idx = build_parent_pool_Aga(fitness_sorted, elite_num, r)
popsize = numel(fitness_sorted);
idx_all = (1:popsize)';
parents = idx_all(1:elite_num);

if elite_num < popsize
    rest_idx = idx_all(elite_num + 1:end);
    mask = rand(numel(rest_idx), 1) < r;
    parents = [parents; rest_idx(mask)];
end

if numel(parents) < 2
    parents = idx_all(1:min(2, popsize));
end

parent_pool_idx = parents;
end

function child = aga_single_point_crossover(p1, p2, N)
M = numel(p1);
if M == 1
    child = p1;
    return;
end
point = randi([1, M - 1]);
child = [p1(1:point), p2(point + 1:end)];
child = unique_fix(child, N);
end

function mutated = aga_point_mutation(ind, N)
M = numel(ind);
mutated = ind;
pos = randi(M);
available = setdiff(1:N, mutated);
if ~isempty(available)
    mutated(pos) = available(randi(numel(available)));
end
mutated = unique_fix(mutated, N);
end

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
coords = wf.candidate_points(ind, :);
cable  = routing_fn(coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, coords);
rt = toc(tic_routing);
end

function fixed = unique_fix(ind, N)
M = numel(ind);
ind = max(1, min(N, ind));
u = unique(ind, 'stable');

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