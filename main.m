close all; clc; clearvars; format long;

%% ========================================================================
% GitHub-ready main entry script
%
% Recommended repository structure:
%
% project_root/
% ├── main.m
% ├── algorithms/
% ├── problems/
% ├── utils/
% ├── Results/     (optional historical results)
% └── results/     (generated experiment outputs)
%% ========================================================================

%% ========================================================================
% Initialize environment
%% ========================================================================
project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end

cd(project_root);
addpath(genpath(project_root));

%% ========================================================================
% Configuration
%% ========================================================================
cfg = struct();

% Experiment parameters
cfg.popsize          = 30;
cfg.max_it           = 100;
cfg.runTime          = 1;
cfg.enable_parallel  = false;

% Deterministic random seed
cfg.base_seed        = 20260316;

% Logging verbosity
cfg.verbose          = true;

% Parallel settings
cfg.parallel_pool_size    = [];
cfg.parallel_idle_timeout = 120;

% Historical results directory
cfg.history_root = fullfile(project_root, 'Results');

% Output directory
cfg.output_root  = fullfile(project_root, 'results');

% Initialization settings
cfg.use_history_pop0 = true;
cfg.allow_history_wf_override = true;
cfg.allow_random_init_when_history_missing = true;

% Save run summaries
cfg.save_run_summary = true;

%% Wind farm cases
cfg.case_list = { ...
    % 'China_Zhuhai_Guishan_Hai', ...
    'Netherlands_Egmond_aan_Zee', ...
    'China_Shanghai_Lingang', ...
    'Netherlands_Prinses_Amaliawindpark', ...
    'Denmark_Nysted', ...
    'UK_Sheringham_Shoal', ...
    'Denmark_Rodsand_II', ...
    'UK_London_Array' ...
};

%% Algorithm registry
cfg.algorithms = { ...
    @GGA, @GA, @AGA, @BPSO, @AGPSO, @BDE, @SaOFGDE, @DOLSSA, @RLPS_TLBO, @EJAYA ...
};

cfg.algonames = { ...
    'GGA','GA','AGA','BPSO','AGPSO','BDE','SaOFGDE','DOLSSA','RLPS_TLBO','EJAYA' ...
};

%% Routing mode
cfg.routing_fn = @cr_sector;   % @cr_mst; @cr_sweep; @cr_sector

%% ========================================================================
% Validate configuration
%% ========================================================================
validate_config(cfg);

%% ========================================================================
% Create output directories
%% ========================================================================
nowStr       = datestr(now, 'yyyymmdd_HHMMSS');
session_name = ['results_' nowStr];

results_dir  = fullfile(cfg.output_root, session_name);
logs_dir     = fullfile(results_dir, 'logs');
meta_dir     = fullfile(results_dir, 'metadata');

safe_mkdir(cfg.output_root);
safe_mkdir(results_dir);
safe_mkdir(logs_dir);
safe_mkdir(meta_dir);

% Save configuration snapshot
config_snapshot = cfg; %#ok<NASGU>
save(fullfile(meta_dir, 'config_snapshot.mat'), 'config_snapshot');
write_text_file(fullfile(meta_dir, 'config_snapshot.txt'), evalc('disp(cfg)'));

%% ========================================================================
% Record environment information
%% ========================================================================
env = collect_environment_info(project_root); %#ok<NASGU>
save(fullfile(meta_dir, 'environment_info.mat'), 'env');
write_text_file(fullfile(meta_dir, 'environment_info.txt'), evalc('disp(env)'));

%% ========================================================================
% Initialize parallel pool if enabled
%% ========================================================================
if cfg.enable_parallel
    ensure_parallel_pool(cfg.parallel_pool_size, cfg.parallel_idle_timeout);
end

%% ========================================================================
% Pre-load problems
%% ========================================================================
nCase = numel(cfg.case_list);

wf_list           = cell(nCase, 1);
turbine_list      = cell(nCase, 1);
wind_name_list    = cell(nCase, 1);
turbine_name_list = cell(nCase, 1);

case_load_status  = false(nCase, 1);
case_load_message = strings(nCase, 1);

for c = 1:nCase
    case_name = cfg.case_list{c};

    try
        [wf, turbine, wind_name, turbine_name] = load_problem_poisson(case_name);

        wf_list{c}           = wf;
        turbine_list{c}      = turbine;
        wind_name_list{c}    = wind_name;
        turbine_name_list{c} = turbine_name;

        case_load_status(c)  = true;
        case_load_message(c) = "OK";

        if cfg.verbose
            fprintf('[INFO] Problem loaded: %s\n', case_name);
        end

    catch ME
        case_load_status(c)  = false;
        case_load_message(c) = string(ME.message);

        if cfg.verbose
            fprintf('[ERROR] Failed to load case: %s\n', case_name);
            fprintf('[ERROR] %s\n', ME.message);
        end
    end
end

problem_summary = table( ...
    string(cfg.case_list(:)), ...
    case_load_status(:), ...
    case_load_message(:), ...
    'VariableNames', {'CaseName', 'Loaded', 'Message'});

writetable(problem_summary, fullfile(meta_dir, 'problem_loading_summary.csv'));

%% ========================================================================
% Main experiment loop
%% ========================================================================
all_run_records = {};

for c = 1:nCase
    case_name = cfg.case_list{c};

    if ~case_load_status(c)
        continue;
    end

    wf      = wf_list{c};
    turbine = turbine_list{c};

    case_result_dir = fullfile(results_dir, sanitize_filename(case_name));
    safe_mkdir(case_result_dir);

    %% Prepare initial populations for all runs
    [Pop0_all, wf, history_info] = prepare_initial_population( ...
        cfg, wf, turbine, case_name, cfg.runTime, cfg.popsize);

    write_text_file( ...
        fullfile(case_result_dir, 'history_init_report.txt'), ...
        history_info.report_text);

    n_turbines   = turbine.turbine_num;
    n_candidates = size(wf.candidate_points, 1);

    %% Algorithm loop
    for a = 1:numel(cfg.algorithms)
        alg      = cfg.algorithms{a};
        algoname = cfg.algonames{a};

        alg_case_dir = fullfile(case_result_dir, sanitize_filename(algoname));
        safe_mkdir(alg_case_dir);

        if cfg.verbose
            fprintf('\n====================================================\n');
            fprintf('Case        : %s\n', case_name);
            fprintf('Algorithm   : %s\n', algoname);
            fprintf('Turbines    : %d\n', n_turbines);
            fprintf('Candidates  : %d\n', n_candidates);
            fprintf('Runs        : %d\n', cfg.runTime);
            fprintf('Iterations  : %d\n', cfg.max_it);
            fprintf('Routing     : %s\n', func2str(cfg.routing_fn));
            fprintf('====================================================\n');
        end

        run_records = cell(cfg.runTime, 1);

        if cfg.enable_parallel
            parfor r = 1:cfg.runTime
                run_records{r} = execute_single_run( ...
                    alg, algoname, case_name, wf, turbine, ...
                    cfg.max_it, r, cfg.popsize, alg_case_dir, Pop0_all{r}, ...
                    cfg.base_seed, cfg.save_run_summary, cfg.routing_fn);
            end
        else
            for r = 1:cfg.runTime
                run_records{r} = execute_single_run( ...
                    alg, algoname, case_name, wf, turbine, ...
                    cfg.max_it, r, cfg.popsize, alg_case_dir, Pop0_all{r}, ...
                    cfg.base_seed, cfg.save_run_summary, cfg.routing_fn);
            end
        end

        run_table = convert_run_records_to_table(run_records);
        writetable(run_table, fullfile(alg_case_dir, 'run_summary.csv'));

        all_run_records = [all_run_records; run_records(:)]; %#ok<AGROW>
    end
end

%% ========================================================================
% Global summary
%% ========================================================================
global_table = convert_run_records_to_table(all_run_records);
writetable(global_table, fullfile(results_dir, 'global_run_summary.csv'));

if cfg.verbose
    fprintf('\nAll experiments finished.\n');
    fprintf('Results stored in: %s\n', results_dir);
end

%% ========================================================================
% Utility functions
%% ========================================================================

function validate_config(cfg)
assert(isstruct(cfg), 'Configuration must be a struct.');
assert(numel(cfg.algorithms) == numel(cfg.algonames), ...
    'Algorithms and names must have the same length.');

for i = 1:numel(cfg.algorithms)
    assert(isa(cfg.algorithms{i}, 'function_handle'), ...
        'Each algorithm entry must be a function handle.');
end

assert(isfield(cfg, 'routing_fn') && isa(cfg.routing_fn, 'function_handle'), ...
    'cfg.routing_fn must be a valid function handle.');
end

function safe_mkdir(folder)
if ~exist(folder, 'dir')
    mkdir(folder);
end
end

function write_text_file(path, txt)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot write file: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', txt);
end

function env = collect_environment_info(project_root)
env = struct();
env.project_root   = project_root;
env.matlab_version = version;
env.computer       = computer;
env.timestamp      = string(datetime);

try
    env.numcores = feature('numcores');
catch
    env.numcores = NaN;
end
end

function name = sanitize_filename(name)
name = char(string(name));
name = regexprep(name, '[<>:"/\\|?*]+', '_');
end

function ensure_parallel_pool(pool_size, idle_timeout)
p = gcp('nocreate');
if ~isempty(p)
    return;
end

pc = parcluster('local');

if isempty(pool_size)
    try
        pool_size = min(pc.NumWorkers, feature('numcores'));
    catch
        pool_size = pc.NumWorkers;
    end
end

try
    pc.IdleTimeout = idle_timeout;
catch
end

parpool(pc, pool_size);
end

function [Pop0_all, wf, history_info] = prepare_initial_population(cfg, wf, turbine, case_name, runTime, popsize)

Pop0_all = cell(runTime, 1);

history_info = struct();
history_info.used_history = false;
history_info.history_file = "";
history_info.report_text  = "";

casePath = fullfile(cfg.history_root, case_name);
matFiles = dir(fullfile(casePath, '*_run*.mat'));

useExisting = false;

if cfg.use_history_pop0 && ~isempty(matFiles)
    try
        fpath = fullfile(matFiles(1).folder, matFiles(1).name);
        S = load(fpath, 'wf', 'generations');

        if cfg.allow_history_wf_override && isfield(S, 'wf')
            wf = S.wf;
        end

        if isfield(S, 'generations') && ~isempty(S.generations) ...
                && isfield(S.generations{1}, 'population')
            refPop = S.generations{1}.population;
            for r = 1:runTime
                Pop0_all{r} = refPop;
            end
            useExisting = true;
            history_info.used_history = true;
            history_info.history_file = string(fpath);
            history_info.report_text = sprintf( ...
                'Case: %s\nInitialization mode: historical population\nHistory file: %s\n', ...
                case_name, fpath);
        end
    catch ME
        history_info.report_text = sprintf( ...
            'Case: %s\nInitialization mode: history loading failed\nMessage: %s\n', ...
            case_name, ME.message);
    end
end

if ~useExisting
    M = turbine.turbine_num;
    N = wf.N_candidate;

    for r = 1:runTime
        rng(cfg.base_seed + r + csum(case_name), 'twister');
        Pop0_all{r} = zeros(popsize, M);
        for i = 1:popsize
            Pop0_all{r}(i, :) = randperm(N, M);
        end
    end

    history_info.used_history = false;
    history_info.report_text = sprintf( ...
        'Case: %s\nInitialization mode: random population\n', case_name);
end
end

function s = csum(str)
s = sum(double(char(str)));
end

function record = execute_single_run( ...
    alg, algoname, case_name, wf, turbine, max_it, run_id, popsize, ...
    output_dir, Pop0, seed, save_summary, routing_fn)

rng(seed + run_id + csum(case_name) + csum(algoname), 'twister');

record = struct();
record.case_name = string(case_name);
record.algorithm = string(algoname);
record.run_id    = run_id;
record.success   = false;
record.best_fit  = NaN;
record.error_message = "";

try
    best_fit = feval(alg, wf, turbine, max_it, run_id, popsize, ...
        algoname, output_dir, Pop0, routing_fn);

    record.success  = true;
    record.best_fit = best_fit;

catch ME
    record.success       = false;
    record.best_fit      = NaN;
    record.error_message = string(getReport(ME, 'basic', 'hyperlinks', 'off'));
end

if save_summary
    save(fullfile(output_dir, sprintf('run_%d_summary.mat', run_id)), 'record');
end
end

function T = convert_run_records_to_table(records)

records = records(~cellfun(@isempty, records));
n = numel(records);

CaseName     = strings(n, 1);
Algorithm    = strings(n, 1);
RunID        = zeros(n, 1);
Success      = false(n, 1);
BestFit      = nan(n, 1);
ErrorMessage = strings(n, 1);

for i = 1:n
    rec = records{i};
    CaseName(i)  = rec.case_name;
    Algorithm(i) = rec.algorithm;
    RunID(i)     = rec.run_id;
    Success(i)   = rec.success;
    BestFit(i)   = rec.best_fit;

    if isfield(rec, 'error_message')
        ErrorMessage(i) = string(rec.error_message);
    else
        ErrorMessage(i) = "";
    end
end

T = table(CaseName, Algorithm, RunID, Success, BestFit, ErrorMessage);
end