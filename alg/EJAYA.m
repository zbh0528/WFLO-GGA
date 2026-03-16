function best_fit = EJAYA(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2021> <single> <integer> <none>
% Enhanced Jaya algorithm
% fi1 --- xxxx --- Random weight controlling attraction toward the best solution
% fi2 --- xxxx --- Random weight controlling repulsion from the worst solution
%
%------------------------------- Reference --------------------------------
% Y. Zhang, A. Chi, and S. Mirjalili. Enhanced Jaya algorithm: A simple but
% efficient optimization method for constrained engineering design
% problems. Knowledge-Based Systems, 2021, 233: 107555.
%
% EJAYA: enhanced JAYA algorithm
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

%% ======================== Basic dimensions ========================
M = turbine.turbine_num;
N = wf.N_candidate;

if nargin < 9 || isempty(routing_fn)
    error('EJAYA:MissingRoutingFunction', ...
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

% Continuous latent representation
X   = double(pop);
old = X;

%% ======================== 2) Initial evaluation ========================
fprintf('Initializing EJAYA population evaluation...\n');
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
fprintf('[EJAYA] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting EJAYA evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    % Current best and worst
    [~, best_idx]  = min(fitness);
    [~, worst_idx] = max(fitness);
    Best  = X(best_idx, :);
    Worst = X(worst_idx, :);
    meanX = mean(X, 1);

    % Weight factors
    fi1 = rand;
    go1 = 1 - fi1;
    fi2 = rand;
    go2 = 1 - fi2;

    % Possible update of old population
    if rand < rand
        old = X;
    end
    old = old(randperm(popsize), :);

    % New generation containers
    new_pop    = pop;
    new_fit    = zeros(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);

    for i = 1:popsize
        ULP = go1 * Best  + fi1 * meanX - X(i, :);
        DLP = go2 * Worst + fi2 * meanX - X(i, :);

        if rand < rand
            Xnew = X(i, :) + rand(1, M) .* ULP - rand(1, M) .* DLP;
        else
            Xnew = X(i, :) + randn(1, M) .* (old(i, :) - X(i, :));
        end

        % Boundary repair
        Xnew = min(max(Xnew, 1), N);

        % Decode into index chromosome
        cand = unique_fix(round(Xnew), N);

        % Evaluate candidate
        [f_new, cf_new, aep_new, cab_new, rt] = ...
            evaluate_individual(cand, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;

        % Acceptance rule
        if f_new < fitness(i)
            X(i, :)        = Xnew;
            pop(i, :)      = cand;
            new_pop(i, :)  = cand;
            new_fit(i)     = f_new;
            new_CF(i)      = cf_new;
            new_AEP(i)     = aep_new;
            new_cables{i}  = cab_new;
        else
            new_pop(i, :)  = pop(i, :);
            new_fit(i)     = fitness(i);
            new_CF(i)      = CF_array(i);
            new_AEP(i)     = AEP_array(i);
            new_cables{i}  = cables{i};
        end
    end

    % Update population state
    pop       = new_pop;
    fitness   = new_fit;
    CF_array  = new_CF;
    AEP_array = new_AEP;
    cables    = new_cables;

    % Save this generation
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
    generations{gen} = generation_struct;

    % Logging
    [best_fit, ~] = min(fitness);
    fprintf('[EJAYA] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('EJAYA finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function fixed = unique_fix(ind, N)
ind = max(min(round(ind), N), 1);
[~, ia] = unique(ind, 'stable');
dup_idx = setdiff(1:numel(ind), ia);
fixed   = ind;

if ~isempty(dup_idx)
    pool = setdiff(1:N, ind(ia));
    need = numel(dup_idx);
    if numel(pool) < need
        extra = setdiff(1:N, ind(ia));
        pool = unique([pool, extra], 'stable');
    end
    fixed(dup_idx) = datasample(pool, need, 'Replace', false);
end
end

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
layout_coords = wf.candidate_points(ind, :);
cable = routing_fn(layout_coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, layout_coords);
rt = toc(tic_routing);
end