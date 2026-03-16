function best_fit = SaOFGDE(wf, turbine, max_it, runtime, popsize, algname, results_dir, Pop0, routing_fn)
% <2025> <single> <integer> <none>
% Self-adaptive opposition-based fractional-order guided differential evolution
% fo_rate --- 0.995 --- Fractional-order memory decay factor in historical guidance
% eps_val --- 0.01  --- Small constant to avoid zero division in success-rate update
% a      --- 0.05  --- Weight controlling contribution of fractional-order historical guidance
% pro1   --- 0.50  --- Initial probability for selecting low-CR mutation strategy
%
%------------------------------- Reference --------------------------------
% Y. Zhang, Z. Zhang, R. Zhong, J. Yu, E. H. Houssein, J. Zhao, and Z. Gao.
% Under complex wind scenarios: Considering large-scale wind turbines in
% wind farm layout optimization via self-adaptive optimal fractional-order
% guided differential evolution. Energy, 2025, 323: 135866.
%
% SaOFGDE: self-adaptive opposition-based fractional-order guided differential evolution
% - Chromosome: a length-M vector of distinct candidate-point indices in 1..N
% - Objective: minimize LCOE
% - Supports a unified initial population Pop0
% - Supports configurable routing function routing_fn
% - Saves full generation history

%% ======================== Basic dimensions ========================
M = turbine.turbine_num;
N = wf.N_candidate;

if nargin < 9 || isempty(routing_fn)
    error('SaOFGDE:MissingRoutingFunction', ...
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
fitness   = zeros(popsize, 1);
CF_array  = zeros(popsize, 1);
AEP_array = zeros(popsize, 1);
cables    = cell(popsize, 1);

fprintf('Initializing SaOFGDE population evaluation...\n');
routing_time_this_gen = 0;
for i = 1:popsize
    [fitness(i), CF_array(i), AEP_array(i), cables{i}, rt] = ...
        evaluate_individual(pop(i, :), wf, turbine, routing_fn);
    fprintf('Initial individual %2d: fitness = %.6f\n', i, fitness(i));
    routing_time_this_gen = routing_time_this_gen + rt;
end

%% ======================== 3) History record ========================
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
fprintf('[SaOFGDE] Generation %3d: best = %.8f | routing = %.2fs\n', ...
    1, best_fit, routing_time_this_gen);

%% ======================== 4) Algorithm parameters ========================
G       = 1;
pNum    = max(1, round(0.1 * popsize));
fo_rate = 0.995;
eps_val = 0.01;
a       = 0.05;

CR   = zeros(popsize, 1);
flag = zeros(popsize, 1);
pro1 = 0.5;
Pop1 = cell(4, 1);

%% ======================== 5) Evolution loop ========================
fprintf('\nStarting SaOFGDE evolution...\n');

for gen = 2:max_it
    gen_tic = tic;
    routing_time_this_gen = 0;

    ns1 = 0; nf1 = 0;
    ns2 = 0; nf2 = 0;

    urnd = rand(popsize, 1);
    F    = 0.1 + 0.9 * rand;

    % --- Adaptive CR ---
    for i = 1:popsize
        if urnd(i) <= pro1
            CR(i)   = 0.05 + 0.1 * rand;
            flag(i) = 1;
        else
            CR(i)   = 0.9 + 0.1 * rand;
            flag(i) = 2;
        end
    end

    % --- New generation containers ---
    new_pop    = zeros(popsize, M);
    new_fit    = zeros(popsize, 1);
    new_CF     = zeros(popsize, 1);
    new_AEP    = zeros(popsize, 1);
    new_cables = cell(popsize, 1);
    V          = zeros(popsize, M);
    U          = zeros(popsize, M);

    % --- Differential evolution ---
    for i = 1:popsize
        [~, indexsortP] = sort(fitness, 'ascend');
        [~, sortCR]     = sort(CR, 'ascend');

        CR_sorted = CR;
        for k = 1:popsize
            CR_sorted(indexsortP(k)) = CR(sortCR(k));
        end
        CR = CR_sorted;

        P = pop(indexsortP, :);

        k0 = randi([1, pNum]);
        P1 = P(k0, :);

        mid_low  = pNum + 1;
        mid_high = popsize - pNum;
        if mid_low > mid_high
            mid_low  = 1;
            mid_high = popsize;
        end
        k1 = randi([mid_low, mid_high]);
        P2 = P(k1, :);

        k2_low = max(1, popsize - pNum + 1);
        k2 = randi([k2_low, popsize]);
        P3 = P(k2, :);

        % --- Historical Pop1 update ---
        if gen <= 4
            Pop1{gen} = P1;
        else
            Pop1{1} = Pop1{2};
            Pop1{2} = Pop1{3};
            Pop1{3} = Pop1{4};
            Pop1{4} = P1;
        end

        % --- Fractional-order guidance ---
        if gen >= ceil((2/3) * max_it) && gen >= 5 && all(~cellfun(@isempty, Pop1))
            P2old = a * ( ...
                (1 / gamma(2) * fo_rate) * Pop1{4} + ...
                (1 / gamma(3) * fo_rate * (1 - fo_rate)) * Pop1{3} + ...
                (1 / gamma(4) * fo_rate * (1 - fo_rate) * (2 - fo_rate)) * Pop1{2} + ...
                (1 / gamma(5) * fo_rate * (1 - fo_rate) * (2 - fo_rate) * (3 - fo_rate)) * Pop1{1});
            P2 = (1 - a) * P2 + P2old;
        end

        V(i, :) = unique_fix(round(P2 + F .* (P1 - P3)), N);

        jrand = randi([1, M]);
        U(i, :) = P(i, :);
        for j = 1:M
            if rand <= CR(i) || j == jrand
                U(i, j) = V(i, j);
            end
        end
        U(i, :) = unique_fix(U(i, :), N);

        new_pop(i, :) = U(i, :);
        [new_fit(i), new_CF(i), new_AEP(i), new_cables{i}, rt] = ...
            evaluate_individual(new_pop(i, :), wf, turbine, routing_fn);
        routing_time_this_gen = routing_time_this_gen + rt;
    end

    % --- Selection ---
    for i = 1:popsize
        if new_fit(i) < fitness(i)
            pop(i, :)    = new_pop(i, :);
            fitness(i)   = new_fit(i);
            CF_array(i)  = new_CF(i);
            AEP_array(i) = new_AEP(i);
            cables{i}    = new_cables{i};

            if flag(i) == 1
                ns1 = ns1 + 1;
            elseif flag(i) == 2
                ns2 = ns2 + 1;
            end
        else
            if flag(i) == 1
                nf1 = nf1 + 1;
            elseif flag(i) == 2
                nf2 = nf2 + 1;
            end
        end
    end

    % --- Update probability ---
    if ns1 + nf1 == 0 || ns2 + nf2 == 0
        pro1 = 0.5;
    else
        s1  = ns1 / (ns1 + nf1) + eps_val;
        s2  = ns2 / (ns2 + nf2) + eps_val;
        ps1 = s1 / (s1 + s2);
        pro1 = ((G - 1) * pro1 + ps1) / G;
    end
    G = G + 1;

    % --- Save this generation ---
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

    % --- Logging ---
    [best_fit, ~] = min(fitness);
    fprintf('[SaOFGDE] Generation %3d: best = %.8f | total = %.2fs | evolution = %.2fs | routing = %.2fs\n', ...
        gen, best_fit, total_time, evolution_time, routing_time_this_gen);
end

%% ======================== 6) Save results ========================
save(fullfile(mat_dir, sprintf('%s_run%02d.mat', algname, runtime)), ...
    'generations', 'wf', 'turbine');

fprintf('SaOFGDE finished. Best LCOE = %.10f\n', best_fit);

end

%% ======================== Helper functions ========================

function [fit, cf, aep, cab, rt] = evaluate_individual(ind, wf, turbine, routing_fn)
tic_routing = tic;
layout_coords = wf.candidate_points(ind, :);
cable = routing_fn(layout_coords, wf);
[fit, cf, aep, cab] = evaluate(wf, turbine, cable, layout_coords);
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