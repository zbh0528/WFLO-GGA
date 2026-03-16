function best_fit = BDE(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2025> <single> <integer> <none>
% Ranking/half-split driven differential evolution
% p_base        --- 0.50 --- Base ratio controlling good/worst population partition
% topK_pBest    --- 10   --- Number of top individuals used as pBest guidance
% use_pBest_mix --- 5    --- Maximum number of pBest candidates used in crossover donor selection
%
%------------------------------- Reference --------------------------------
% J. Li, Z. Zhang, T. Zheng, J. Tang, Z. Lei, and S. Gao. Discrete
% bi-population differential evolution for optimizing complex wind farm
% layouts in diverse terrains. Energy, 2025: 137885.
%
% BDE: ranking/half-split driven differential evolution
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

%% ======================== Parameter preparation ========================
M = turbine.turbine_num;
N = wf.N_candidate;

% BDE parameters
p_base        = 0.50;
topK_pBest    = 10;
use_pBest_mix = 5;

if nargin < 9 || isempty(routing_fn)
    error('BDE:MissingRoutingFunction', ...
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
fprintf('Initializing BDE population evaluation...\n');
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
fprintf('[BDE] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting BDE evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    % --- Adaptive parameters ---
    pjd     = gen / max_it;
    self_sy = 2^exp(1 - max_it / (max_it + 1 - gen));
    F       = 0.5 * self_sy;
    CR      = 0.1 * pjd;
    p       = p_base;

    % --- Sort population by fitness (ascending is better) ---
    [fitness, order] = sort(fitness, 'ascend');
    pop       = pop(order, :);
    CF_array  = CF_array(order);
    AEP_array = AEP_array(order);
    cables    = cables(order);

    % --- Partition population ---
    N_half  = round(popsize / 2);
    K_pBest = min(topK_pBest, popsize);
    pBest   = pop(1:K_pBest, :);

    goodPops  = pop(1:N_half, :);
    worstPops = pop(N_half+1:end, :);

    gP = goodPops(randperm(N_half), :);
    wP = worstPops(randperm(size(worstPops, 1)), :);

    idx_cut = max(0, min(N_half, floor(p * N_half)));
    gp1 = gP(1:idx_cut, :);
    gp2 = gP(idx_cut+1:N_half, :);

    if size(wP, 1) < N_half
        wP = repmat(wP, ceil(N_half / size(wP, 1)), 1);
    end
    wP  = wP(1:N_half, :);
    wp1 = wP(1:N_half-idx_cut, :);
    wp2 = wP(N_half-idx_cut+1:N_half, :);

    V_pops = [gp1; wp1];
    U_pool = [gp2; wp2];

    % --- Differential mutation ---
    r1 = randperm(N_half);
    r2 = randperm(N_half);
    r3 = randperm(N_half);
    V  = V_pops(r1, :) + F * (V_pops(r2, :) - U_pool(r3, :));

    % --- Binomial crossover ---
    trial = zeros(N_half, M);
    for i = 1:N_half
        j_rand = randi(M);
        mask = rand(1, M) < CR;
        mask(j_rand) = true;

        pBest_pick_max = min(use_pBest_mix, size(pBest, 1));
        donor = pBest(randi(pBest_pick_max), :);

        trial_real = mask .* V(i, :) + (~mask) .* donor;
        trial(i, :) = discretize_and_fix(trial_real, N, M);
    end

    % --- Evaluate trial individuals ---
    trial_fit  = zeros(N_half, 1);
    trial_CF   = zeros(N_half, 1);
    trial_AEP  = zeros(N_half, 1);
    trial_cabs = cell(N_half, 1);

    for i = 1:N_half
        [trial_fit(i), trial_CF(i), trial_AEP(i), trial_cabs{i}, rt] = ...
            evaluate_individual(trial(i, :), wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;
    end

    % --- Selection and replacement ---
    for i = 1:N_half
        if trial_fit(i) < fitness(i)
            pop(i, :)    = trial(i, :);
            fitness(i)   = trial_fit(i);
            CF_array(i)  = trial_CF(i);
            AEP_array(i) = trial_AEP(i);
            cables{i}    = trial_cabs{i};
        end
    end

    % --- Record ---
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

    % --- Logging ---
    [best_fit, ~] = min(fitness);
    fprintf('[BDE] Generation %3d: best = %.8f | F = %.3f | CR = %.3f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, F, CR, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');
fprintf('BDE finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function idx_vec = discretize_and_fix(real_vec, N, M)
cand = round(real_vec);
cand = max(1, min(N, cand));
idx_vec = unique_fix(cand, N);

if numel(idx_vec) > M
    idx_vec = idx_vec(1:M);
elseif numel(idx_vec) < M
    pool = setdiff(1:N, idx_vec);
    addn = M - numel(idx_vec);
    if numel(pool) < addn
        pool = [pool, setdiff(1:N, idx_vec)];
        pool = unique(pool, 'stable');
    end
    add = datasample(pool, addn, 'Replace', false);
    idx_vec = [idx_vec, add];
end
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

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
coords = wf.candidate_points(ind, :);
cable  = routing_fn(coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, coords);
rt = toc(tic_routing);
end