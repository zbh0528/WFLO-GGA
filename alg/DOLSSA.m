function best_fit = DOLSSA(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2024> <single> <integer> <none>
% Dynamic opposition-based learning enhanced sparrow search algorithm
% DiscovererRatio --- 0.20 --- Proportion of discoverers in the population
% CompanionRatio  --- 0.70 --- Proportion of companions in the population
%
%------------------------------- Reference --------------------------------
% Y. Zhu, Y. Guo, T. Hu, C. Wu, and L. Zhang. Wind farm layout optimization
% based on dynamic opposite learning‐enhanced sparrow search algorithm.
% International Journal of Energy Research, 2024, 2024(1): 4322211.
%
% DOLSSA: dynamic opposition-based learning sparrow search algorithm
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

M = turbine.turbine_num;
N = wf.N_candidate;

if nargin < 9 || isempty(routing_fn)
    error('DOLSSA:MissingRoutingFunction', ...
        'routing_fn must be provided, e.g., @cr_mst, @cr_sweep, or @cr_sector.');
end

%% ======================== Output folder ========================
mat_dir = results_dir;
if ~exist(mat_dir, 'dir')
    mkdir(mat_dir);
end

%% ======================== 1) Initialization ========================
if nargin >= 8 && ~isempty(Pop0)
    pop = Pop0;
else
    pop = zeros(popsize, M);
    for i = 1:popsize
        pop(i, :) = randperm(N, M);
    end
end

% Generate opposite population
opp_pop = OppositePopulation(pop, wf, 1);
combined_pop = [pop; opp_pop];

% Evaluate combined population
fitness   = zeros(size(combined_pop, 1), 1);
CF_array  = zeros(size(combined_pop, 1), 1);
AEP_array = zeros(size(combined_pop, 1), 1);
cables    = cell(size(combined_pop, 1), 1);
routing_time_this_gen = 0;

for i = 1:size(combined_pop, 1)
    [fitness(i), CF_array(i), AEP_array(i), cables{i}, rt] = ...
        evaluate_individual(combined_pop(i, :), wf, turbine, routing_fn);
    routing_time_this_gen = routing_time_this_gen + rt;
end

% Select the best popsize individuals
[~, sortIdx] = sort(fitness, 'ascend');
pop       = combined_pop(sortIdx(1:popsize), :);
fitness   = fitness(sortIdx(1:popsize));
CF_array  = CF_array(sortIdx(1:popsize));
AEP_array = AEP_array(sortIdx(1:popsize));
cables    = cables(sortIdx(1:popsize));

% Global best
[best_fit, best_idx] = min(fitness);
gBest = struct( ...
    'Position',  pop(best_idx, :), ...
    'Cost',      best_fit, ...
    'CF',        CF_array(best_idx), ...
    'AEP',       AEP_array(best_idx), ...
    'CableInfo', cables{best_idx});

%% ======================== 2) History initialization ========================
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

fprintf('[%s] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    algname, 1, gBest.Cost, routing_time_this_gen);

%% ======================== 3) Main loop ========================
for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    Offdecs = zeros(popsize, M);
    D = M;
    Objs = fitness;

    % Discoverers
    DiscovererSize = max(1, round(0.2 * popsize));
    for i = 1:DiscovererSize
        if rand > 0.8
            Offdecs(i, :) = pop(i, :) .* exp(-i / (rand * max_it + eps));
        else
            Offdecs(i, :) = pop(i, :) + randn(1, D);
        end
    end

    % Companions
    [~, idx_best]  = min(Objs);
    [~, idx_worst] = max(Objs);
    Xbest  = pop(idx_best, :);
    Xworst = pop(idx_worst, :);

    CompanionSize = round(0.7 * popsize);
    companion_end = min(popsize, DiscovererSize + CompanionSize);

    for i = DiscovererSize + 1:companion_end
        if i > popsize / 2
            Offdecs(i, :) = randn(1, D) .* exp((Xworst - pop(i, :)) / (i^2 + eps));
        else
            Offdecs(i, :) = Xbest + abs(pop(i, :) - Xbest) .* (randi([0, 1], 1, D) * 2 - 1);
        end
    end

    % Warners
    for i = companion_end + 1:popsize
        if Objs(i) > min(Objs)
            Offdecs(i, :) = gBest.Position + randn(1, D) .* (pop(i, :) - gBest.Position);
        else
            K = 2 * rand - 1;
            Offdecs(i, :) = pop(i, :) + K * (pop(i, :) - Xworst) / (Objs(i) - max(Objs) + eps);
        end
    end

    % Repair to valid index chromosomes
    for i = 1:popsize
        Offdecs(i, :) = unique_fix(round(Offdecs(i, :)), N);
    end

    % Evaluate offspring
    new_fit    = zeros(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);

    for i = 1:popsize
        [new_fit(i), new_CF(i), new_AEP(i), new_cables{i}, rt] = ...
            evaluate_individual(Offdecs(i, :), wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;
    end

    % Greedy selection
    improved = new_fit < fitness;
    pop(improved, :)     = Offdecs(improved, :);
    fitness(improved)    = new_fit(improved);
    CF_array(improved)   = new_CF(improved);
    AEP_array(improved)  = new_AEP(improved);
    cables(improved)     = new_cables(improved);

    % Update global best
    [cur_best, best_idx] = min(fitness);
    if cur_best < gBest.Cost
        gBest.Cost      = cur_best;
        gBest.Position  = pop(best_idx, :);
        gBest.CF        = CF_array(best_idx);
        gBest.AEP       = AEP_array(best_idx);
        gBest.CableInfo = cables{best_idx};
    end

    % Save history
    total_time = toc(gen_tic);
    generation_struct = struct();
    generation_struct.population        = pop;
    generation_struct.fitness           = fitness;
    generation_struct.CF                = CF_array;
    generation_struct.AEP               = AEP_array;
    generation_struct.cables            = cables;
    generation_struct.runtime           = total_time;
    generation_struct.runtime_routing   = routing_time_this_gen;
    generation_struct.runtime_evolution = total_time - routing_time_this_gen;
    generations{gen} = generation_struct;

    fprintf('[%s] Generation %3d: best = %.8f | total = %.2fs | routing = %.2fs\n', ...
        algname, gen, gBest.Cost, total_time, routing_time_this_gen);

    best_fit = gBest.Cost;
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('[%s] Finished. Best LCOE = %.8f\n', algname, best_fit);
end

%% ======================== Helper functions ========================
function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
coords = wf.candidate_points(ind, :);
cable  = routing_fn(coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, coords);
rt = toc(tic_routing);
end

function OppPop = OppositePopulation(Pop, wf, Case)
[popsize, M] = size(Pop);
N = wf.N_candidate;
OppPop = zeros(popsize, M);

for i = 1:popsize
    opp_ind = zeros(1, M);
    for j = 1:M
        gene = Pop(i, j);
        opp_gene = N - gene + 1;
        if Case == 1
            opp_gene = opp_gene + randn * 0.05 * N;
        else
            opp_gene = opp_gene + randn * 0.01 * N;
        end
        opp_gene = max(min(round(opp_gene), N), 1);
        opp_ind(j) = opp_gene;
    end
    opp_ind = unique_fix(opp_ind, N);
    OppPop(i, :) = opp_ind;
end
end

function fixed = unique_fix(ind, N)
M = numel(ind);
ind = max(min(round(ind), N), 1);
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

add = randsample(pool, need, false);
fixed = [u, add];
end