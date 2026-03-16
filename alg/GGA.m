function best_fit = GGA(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <xxxx> <single> <integer> <none>
% Sector-based adaptive genetic algorithm
% pm --- 0.50 --- Mutation probability
%
%------------------------------- Reference --------------------------------
% xxxx
%
% SAGA: Sector-based Adaptive Genetic Algorithm
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history
%% ======================== Parameter preparation ========================
M = turbine.turbine_num;
N = wf.N_candidate;
candidate_coords = wf.candidate_points;

pm = 0.5;
N_all = 1:N;

if nargin < 9 || isempty(routing_fn)
    error('SAGA:MissingRoutingFunction', ...
        'routing_fn must be provided, e.g., @cr_mst, @cr_sweep, or @cr_sector.');
end

% Candidate-point polar angles in [0, 2*pi)
all_angles = atan2(candidate_coords(:,2) - wf.substation(2), ...
                   candidate_coords(:,1) - wf.substation(1));
all_angles(all_angles < 0) = all_angles(all_angles < 0) + 2 * pi;

% Individual angle parameters
theta0 = rand(popsize, 1) * pi;

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
fprintf('Initializing SAGA population evaluation...\n');
fitness   = zeros(popsize, 1);
CF_array  = zeros(popsize, 1);
AEP_array = zeros(popsize, 1);
cables    = cell(popsize, 1);

routing_time_this_gen = 0;
for i = 1:popsize
    [fitness(i), CF_array(i), AEP_array(i), cables{i}, routing_time] = ...
        evaluate_individual(pop(i, :), wf, turbine, routing_fn);
    fprintf('Initial individual %2d: fitness = %.6f\n', i, fitness(i));
    routing_time_this_gen = routing_time_this_gen + routing_time;
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
generation_struct.diversity         = compute_diversity(pop);
generation_struct.theta0            = theta0;
generations{1} = generation_struct;

[best_fit, best_idx] = min(fitness);
fprintf('[SAGA] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting SAGA evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;
    diversity = compute_diversity(pop);

    new_pop    = zeros(popsize, M);
    new_fit    = inf(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);
    new_theta0 = zeros(popsize, 1);

    for i = 1:popsize
        [p1, p2] = select_parents(pop, fitness);

        theta_i = rand * pi;
        new_theta0(i) = theta_i;

        theta2 = theta_i + pi;
        regionA = get_region(theta_i, theta2, all_angles);
        regionB = setdiff(N_all, regionA);

        [child1, child2] = crossover(p1, p2, regionA, regionB, N);
        child = child1;
        if rand > 0.5
            child = child2;
        end
        child = fix_child_size(unique(child, 'stable'), M, N);

        if rand < pm
            child = mutation(child, N, M);
        end
        child = unique_fix(child, N);

        [fit, cf, aep, cable, routing_time1] = ...
            evaluate_individual(child, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + routing_time1;

        new_pop(i, :)   = child;
        new_fit(i)      = fit;
        new_CF(i)       = cf;
        new_AEP(i)      = aep;
        new_cables{i}   = cable;
    end

    improved = new_fit < fitness;
    pop(improved, :)      = new_pop(improved, :);
    fitness(improved)     = new_fit(improved);
    CF_array(improved)    = new_CF(improved);
    AEP_array(improved)   = new_AEP(improved);
    theta0(improved)      = new_theta0(improved);
    cables(improved)      = new_cables(improved);

    [best_fit, best_idx] = min(fitness);

    total_time     = toc(gen_tic);
    evolution_time = total_time - routing_time_this_gen;

    generation_struct = struct();
    generation_struct.population        = pop;
    generation_struct.fitness           = fitness;
    generation_struct.CF                = CF_array;
    generation_struct.AEP               = AEP_array;
    generation_struct.cables            = cables;
    generation_struct.runtime           = total_time;
    generation_struct.runtime_routing   = routing_time_this_gen;
    generation_struct.runtime_evolution = evolution_time;
    generation_struct.diversity         = diversity;
    generation_struct.theta0            = theta0;
    generations{gen} = generation_struct;

    fprintf('[SAGA] Generation %3d: best = %.8f | diversity = %.4f | total = %.2fs | evolution = %.2fs | routing = %.2fs | theta0 = %.1f deg\n', ...
        gen, best_fit, diversity, total_time, evolution_time, routing_time_this_gen, rad2deg(theta0(best_idx)));
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');
fprintf('SAGA finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function [p1, p2] = select_parents(pop, fitness)
idx1 = roulette_selection(fitness);
idx2 = roulette_selection(fitness);
if idx1 == idx2
    idx2 = mod(idx2 + randi([1, size(pop, 1) - 1]), size(pop, 1));
    if idx2 == 0
        idx2 = size(pop, 1);
    end
end
p1 = pop(idx1, :);
p2 = pop(idx2, :);
end

function region = get_region(theta1, theta2, all_angles)
if theta1 < theta2
    region = find(all_angles >= theta1 & all_angles < theta2);
else
    region = find(all_angles >= theta1 | all_angles < theta2);
end
end

function [child1, child2] = crossover(p1, p2, regionA, regionB, N)
maskA = false(1, N);
maskB = false(1, N);
maskA(regionA) = true;
maskB(regionB) = true;
child1 = [p1(maskA(p1)), p2(maskB(p2))];
child2 = [p2(maskA(p2)), p1(maskB(p1))];
child1 = unique(child1, 'stable');
child2 = unique(child2, 'stable');
end

function fixed = fix_child_size(child, M, N)
child = unique(child, 'stable');
len = numel(child);

if len < M
    mark = false(1, N);
    mark(child) = true;
    available = find(~mark);
    if numel(available) >= M - len
        fixed = [child, available(randperm(numel(available), M - len))];
    else
        fixed = [child, available];
        fixed = unique_fix(fixed, N);
        if numel(fixed) > M
            fixed = fixed(1:M);
        end
    end
elseif len > M
    fixed = child(randperm(len, M));
else
    fixed = child;
end
fixed = unique_fix(fixed, N);
end

function child = mutation(child, N, M)
idx = randi(M);
mark = false(1, N);
mark(child) = true;
available = find(~mark);
if ~isempty(available)
    child(idx) = available(randi(numel(available)));
end
child = unique_fix(child, N);
end

function div = compute_diversity(pop)
[popsize, ~] = size(pop);
N = max(pop(:));
bin = false(popsize, N);
for i = 1:popsize
    bin(i, pop(i, :)) = true;
end
div_sum = 0;
for i = 1:popsize-1
    A = bin(i, :);
    for j = i+1:popsize
        B = bin(j, :);
        inter = sum(A & B);
        union = sum(A | B);
        if union == 0
            diff = 0;
        else
            diff = 1 - inter / union;
        end
        div_sum = div_sum + diff;
    end
end
pair_count = popsize * (popsize - 1) / 2;
if pair_count == 0
    div = 0;
else
    div = div_sum / pair_count;
end
end

function idx = roulette_selection(fitness)
epsilon = 1e-8;
scores = max(fitness) - fitness + epsilon;
cdf = cumsum(scores);
r = rand * cdf(end);
idx = find(cdf >= r, 1);
end

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
coords = wf.candidate_points(ind, :);
cable = routing_fn(coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, coords);
rt = toc(tic_routing);
end

function fixed = unique_fix(ind, N)
M = numel(ind);
ind = max(1, min(N, round(ind)));
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