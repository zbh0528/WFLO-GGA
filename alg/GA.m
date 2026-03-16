function best_fit = GA(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
<1975> <single> <integer> <none>
% Classical genetic algorithm
% mutation_rate --- 0.05 --- Probability of mutating each gene
% elite_num     --- 1    --- Number of elite individuals preserved each generation
%
%------------------------------- Reference --------------------------------
% J. H. Holland. Adaptation in Natural and Artificial Systems. University
% of Michigan Press, 1975.
%
% GA: classical genetic algorithm
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history
%
% Inputs:
%   wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn
%
% Output:
%   best_fit    best LCOE found in this run

%% ======================== Parameter preparation ========================
M = turbine.turbine_num;
N = wf.N_candidate;

mutation_rate = 0.05;
elite_num     = 1;

if nargin < 9 || isempty(routing_fn)
    error('GA:MissingRoutingFunction', ...
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
fprintf('Initializing GA population evaluation...\n');

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
generations{1} = generation_struct;

[best_fit, ~] = min(fitness);
fprintf('[GA] Generation %3d: best = %.8f | routing time = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 3) Evolution loop ========================
fprintf('\nStarting GA evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    %% Elite preservation
    [~, elite_idx] = min(fitness);
    elite = pop(elite_idx, :);

    elite_data = struct( ...
        'fit',   fitness(elite_idx), ...
        'cf',    CF_array(elite_idx), ...
        'aep',   AEP_array(elite_idx), ...
        'cable', cables{elite_idx});

    %% New population containers
    new_pop    = zeros(popsize, M);
    new_fit    = zeros(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);

    %% Put elite into the first slot
    new_pop(1, :) = elite;
    new_fit(1)    = elite_data.fit;
    new_CF(1)     = elite_data.cf;
    new_AEP(1)    = elite_data.aep;
    new_cables{1} = elite_data.cable;

    %% Selection + crossover + mutation
    idx_new = elite_num + 1;

    while idx_new <= popsize
        p1 = tournament_selection(fitness, 2);
        p2 = tournament_selection(fitness, 2);

        c1 = single_point_crossover(pop(p1, :), pop(p2, :), N);
        c2 = single_point_crossover(pop(p2, :), pop(p1, :), N);

        c1 = simple_point_mutation(c1, N, mutation_rate);
        c2 = simple_point_mutation(c2, N, mutation_rate);

        children = {c1, c2};

        for j = 1:2
            if idx_new > popsize
                break;
            end

            [new_fit(idx_new), new_CF(idx_new), new_AEP(idx_new), ...
                new_cables{idx_new}, routing_time] = ...
                evaluate_individual(children{j}, wf, turbine, routing_fn);

            routing_time_this_gen = routing_time_this_gen + routing_time;
            new_pop(idx_new, :) = children{j};
            idx_new = idx_new + 1;
        end
    end

    %% Replace old population
    pop       = new_pop;
    fitness   = new_fit;
    CF_array  = new_CF;
    AEP_array = new_AEP;
    cables    = new_cables;

    %% Record this generation
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

    %% Logging
    [best_fit, ~] = min(fitness);
    fprintf('[GA] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 4) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('GA finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
layout_coords = wf.candidate_points(ind, :);
cable = routing_fn(layout_coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, layout_coords);
rt = toc(tic_routing);
end

function idx = tournament_selection(fitness, k)
candidates = randperm(length(fitness), k);
[~, best_idx] = min(fitness(candidates));
idx = candidates(best_idx);
end

function child = single_point_crossover(p1, p2, N)
M = length(p1);
point = randi([1, M - 1]);
child = [p1(1:point), p2(point+1:end)];
child = unique_fix(child, N);
end

function mutated = simple_point_mutation(ind, N, rate)
mutated = ind;
for i = 1:length(ind)
    if rand < rate
        available = setdiff(1:N, mutated);
        if ~isempty(available)
            mutated(i) = datasample(available, 1);
        end
    end
end
mutated = unique_fix(mutated, N);
end