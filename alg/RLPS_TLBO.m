function best_fit = RLPS_TLBO(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2024> <single> <integer> <none>
% Reinforcement learning based phase-selection teaching–learning-based optimization
% alpha   --- 0.10 --- Learning rate for Q-table update
% gamma_q --- 0.90 --- Discount factor in reinforcement learning update
% f_values --- [-0.05,0,0.05] --- Action set controlling exploration probability adjustment
%
%------------------------------- Reference --------------------------------
% X. Yu and W. Zhang. A teaching-learning-based optimization algorithm with
% reinforcement learning to address wind farm layout optimization problem.
% Applied Soft Computing, 2024, 151: 111135.
%
% RLPS_TLBO: reinforcement learning based phase-selection TLBO
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

M = turbine.turbine_num;
N = wf.N_candidate;

if nargin < 9 || isempty(routing_fn)
    error('RLPS_TLBO:MissingRoutingFunction', ...
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

%% ======================== 2) Individual evaluation ========================
fitness   = zeros(popsize, 1);
CF_array  = zeros(popsize, 1);
AEP_array = zeros(popsize, 1);
cables    = cell(popsize, 1);
routing_time_this_gen = 0;

for i = 1:popsize
    [fitness(i), CF_array(i), AEP_array(i), cables{i}, rt] = ...
        evaluate_individual(pop(i, :), wf, turbine, routing_fn);
    routing_time_this_gen = routing_time_this_gen + rt;
end

[best_fit, best_idx] = min(fitness);
gBest = struct( ...
    'Position',  pop(best_idx, :), ...
    'Cost',      best_fit, ...
    'CF',        CF_array(best_idx), ...
    'AEP',       AEP_array(best_idx), ...
    'CableInfo', cables{best_idx});

%% ======================== 3) Q-table initialization ========================
Q_table = zeros(2, 3);
alpha   = 0.1;
gamma_q = 0.9;
f_values = [-0.05, 0, 0.05];
f = f_values(randi(3));
F = rand + f;
F = max(0, min(1, F));

%% ======================== 4) History container ========================
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
    algname, 1, best_fit, routing_time_this_gen);

%% ======================== 5) Main loop ========================
for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    for i = 1:popsize
        if rand < F
            % Teacher phase
            TF = randi([1, 2]);
            [~, bestIdx] = min(fitness);
            OffDec = pop(i, :) + rand * (pop(bestIdx, :) - TF * mean(pop, 1));
            state = 1;
        else
            % Learner phase
            r1 = randi(popsize);
            while r1 == i
                r1 = randi(popsize);
            end

            r2 = randi(popsize);
            while r2 == i || r2 == r1
                r2 = randi(popsize);
            end

            Learner = pop(i, :) + rand * (pop(r1, :) - pop(r2, :));
            if rand <= 0.5
                OffDec = Learner;
            else
                OffDec = pop(i, :);
            end
            state = 2;
        end

        OffDec = unique_fix(round(OffDec), N);

        [f_new, cf_new, aep_new, cab_new, rt] = ...
            evaluate_individual(OffDec, wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;

        if f_new < fitness(i)
            pop(i, :)    = OffDec;
            fitness(i)   = f_new;
            CF_array(i)  = cf_new;
            AEP_array(i) = aep_new;
            cables{i}    = cab_new;
            reward = 1;
        else
            reward = 0;
        end

        % Q-table update
        aj = find(abs(f_values - f) < 1e-12, 1);
        if isempty(aj)
            [~, aj] = min(abs(f_values - f));
        end

        Q_table(state, aj) = Q_table(state, aj) + ...
            alpha * (reward + gamma_q * max(Q_table(state, :)) - Q_table(state, aj));

        f = calf(Q_table, state, 1);
        F = F + f;
        F = max(0, min(1, F));
    end

    % Update global best
    [cur_best, best_idx] = min(fitness);
    if cur_best < gBest.Cost
        gBest.Cost      = cur_best;
        gBest.Position  = pop(best_idx, :);
        gBest.CF        = CF_array(best_idx);
        gBest.AEP       = AEP_array(best_idx);
        gBest.CableInfo = cables{best_idx};
    end

    % Save current generation
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

%% ======================== 6) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('[%s] Finished. Best LCOE = %.8f\n', algname, best_fit);
end

%% ======================== Helper functions ========================
function f = calf(Q, state, T)
f_values = [-0.05, 0, 0.05];
Q_values = Q(state, :);

if all(abs(Q_values - Q_values(1)) < 1e-12)
    prob = ones(1, length(Q_values)) / length(Q_values);
else
    exp_values = exp(Q_values / T);
    prob = exp_values / sum(exp_values);
end

chosen_action = randsample(1:length(f_values), 1, true, prob);
f = f_values(chosen_action);
end

function fixed = unique_fix(ind, N)
M = numel(ind);
ind = max(min(ind, N), 1);
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

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
coords = wf.candidate_points(ind, :);
cable = routing_fn(coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, coords);
rt = toc(tic_routing);
end