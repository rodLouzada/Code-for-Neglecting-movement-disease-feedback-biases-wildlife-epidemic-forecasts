function Simulation_Framework()

rng(42, 'multFibonacci');

sir_mode = 1;  % 0 = SI, 1 = SIR

output_folder = 'path\to\folder\MultiLambda Z Score SIR Status Cohesion';

    % World simulation parameters
    board_size_x = 50;  % Board size X axis
    board_size_y = 50;  % Board size Y axis
    simulation_time = 1;  % Simulation time
    time_step = 1;  % Time step
    num_steps = simulation_time / time_step;  % Number of total steps
    num_infected_agents = 1;  % Number of infected agents per group
    num_infected_groups = 1;  % Number of groups with initially infected agents
    pathogen_deposition = 1;  % Pathogen load shed at each time step by infected agents
    num_repeats = 1;  % Number of repetitions to reduce the impact of stochasticity

    num_agents = 1;  % Number of agents per group
    num_groups = 1;   % Number of groups 1, 5, 

    % Disease parameters
    pathogen_decay = 0.25;
    recovery_rate = 0.005;  % Recovery rate
    infection_rate = 0.1;  % Rate of infection

    % Movement parameters
    % Movement function parameter for diffusion   
    diffusion = 3; %1, 2, 3,4,5, 6,7,8,9,10,15,20
    drift = 0.2;  % Movement function parameter for drift
    % LAMBDA Infected movement impact

    infected_movement_impact_c = 0.1; % 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1
    infected_movement_impact_d = 0.1;

    % NEW: Group cohesion coefficients to sweep
    group_cohesion = 0.9; % 0, 0.1, 0.5, 0.9, 0.05, 0.25, 0.75, 0.95

fprintf('Starting Simulation\n');

simulation( ...
           num_agents, ...
           num_group, ...
           board_size_x, ...
           board_size_y, ...
           num_infected_agents, ...
           num_infected_groups, ...
           recovery_rate, ...
           pathogen_deca, ...
           pathogen_deposition, ...
           num_steps, ...
           diffusion, ...
           drift, ...
           time_step, ...
           num_repeats, ...
           infection_rate, ...
           infected_movement_impact_c, ... 
           infected_movement_impact_d,  ... 
           group_cohesion, ...                              
           sir_mode, ...
           output_folder);

                                 
                                   
end


function simulation(num_agents, num_groups, board_size_x, board_size_y, num_infected_agents, num_infected_groups, ...
    recovery_rate, pathogen_decay_rate, pathogen_deposition, num_steps, diffusion, drift, time_step, ...
    num_repeats, infection_rate, lambda_c, lambda_d, group_cohesion, sir_mode, output_folder)

num_per_group = num_agents / num_groups;
burn_in_steps = 10;

txt_SIR       = zeros(num_steps, 3);
txt_Peak      = zeros(num_repeats, 2);
txt_std_SIR   = zeros(num_steps, 3);
txt_prevalence= zeros(num_repeats, 3);

% SIR status tensor: agents x repeats x timesteps+1
SIR_status = zeros(num_agents, num_repeats, num_steps+1);

% Optional coordinate record for rp==1
coord_record = zeros(num_agents * num_steps, 5); % step, ID, X, Y, SIR

dWx = create_random_3d_array(num_agents, num_repeats, num_steps, burn_in_steps);
dWy = create_random_3d_array(num_agents, num_repeats, num_steps, burn_in_steps);

board = create_grid_world_array(board_size_x, board_size_y, num_repeats);

x_init_coord = randomize_group_start(board_size_x, num_agents, num_groups, num_repeats, num_per_group);
y_init_coord = randomize_group_start(board_size_y, num_agents, num_groups, num_repeats, num_per_group);

x_current_coord = x_init_coord;
y_current_coord = y_init_coord;

% Burn-in (no cohesion in burn-in)
for bis = 1:burn_in_steps
    x_current_coord = Next_Sim_Step(diffusion, drift, x_init_coord, x_current_coord, dWx, board_size_x, bis);
    y_current_coord = Next_Sim_Step(diffusion, drift, y_init_coord, y_current_coord, dWy, board_size_y, bis);
end

% Infect initial agents
SIR_status(:,:,1) = infect_agents(num_infected_groups, num_infected_agents, num_per_group, num_agents, num_repeats);

% --- Per-timestep loop ---
for t = 1:num_steps

    % Aggregate S/I/R means
    txt_SIR(t,1) = sum(sum(SIR_status(:,:,t) == 0)) / num_repeats;
    txt_SIR(t,2) = sum(sum(SIR_status(:,:,t) == 1)) / num_repeats;
    txt_SIR(t,3) = sum(sum(SIR_status(:,:,t) == 2)) / num_repeats;

    % Std across repeats (counts per repeat at time t)
    txt_std_SIR(t,1) = std(sum(SIR_status(:,:,t) == 0, 1));
    txt_std_SIR(t,2) = std(sum(SIR_status(:,:,t) == 1, 1));
    txt_std_SIR(t,3) = std(sum(SIR_status(:,:,t) == 2, 1));

    if t == num_steps
        % Prevalence at the last step per repetition
        txt_prevalence(:,1) = sum(SIR_status(:,:,t) == 0, 1); % S
        txt_prevalence(:,2) = sum(SIR_status(:,:,t) == 1, 1); % I
        txt_prevalence(:,3) = sum(SIR_status(:,:,t) == 2, 1); % R
    end

    % ----------- NEW: compute group centers (for current positions) -----------
    % group_center_x/y: (num_groups x num_repeats), hard-edge map
    group_center_x = zeros(num_groups, num_repeats);
    group_center_y = zeros(num_groups, num_repeats);

    for rp = 1:num_repeats
        for g = 1:num_groups
            g_start = (g-1)*num_per_group + 1;
            g_end   = g*num_per_group;

            % mean over the agents belonging to group g in repeat rp
            group_center_x(g, rp) = mean(x_current_coord(g_start:g_end, rp));
            group_center_y(g, rp) = mean(y_current_coord(g_start:g_end, rp));
        end
    end
    % --------------------------------------------------------------------------

    % Track peak (note: replicates original logic)
    for rp = 1:num_repeats
        % FIX
curr_I = sum(SIR_status(:,rp,t) == 1);
if curr_I > txt_Peak(rp,1)
    txt_Peak(rp,1) = curr_I;
    txt_Peak(rp,2) = t;
end

    end

    % Agent updates
    for rp = 1:num_repeats
        for ag = 1:num_agents

            x = round(x_current_coord(ag, rp));
            y = round(y_current_coord(ag, rp));

            % Save coords for rp==1
            if rp == 1
                coord_record((ag-1)*num_steps + t, :) = [t, ag, x_current_coord(ag,rp), y_current_coord(ag,rp), SIR_status(ag,rp,t)];
            end

            % Determine agent's group and its center (for this repeat)
            g = ceil(ag / num_per_group);
            cx = group_center_x(g, rp);
            cy = group_center_y(g, rp);

            if SIR_status(ag,rp,t) == 1
                % Infected: deposit pathogen time-stamp at cell
                board{x,y,rp} = [board{x,y,rp}, t];

                % Recovery check
                if rand() < recovery_rate
                    SIR_status(ag,rp,t+1) = 2;
                else
                    SIR_status(ag,rp,t+1) = 1;
                end

                % Move with λ-modified drift/diffusion + cohesion
                x_current_coord(ag,rp) = Next_Step_Cohesion( ...
                    diffusion*lambda_d, drift*lambda_c, group_cohesion, ...
                    x_init_coord(ag,rp), x_current_coord(ag,rp), cx, ...
                    dWx(ag,rp,t+burn_in_steps), board_size_x);

                y_current_coord(ag,rp) = Next_Step_Cohesion( ...
                    diffusion*lambda_d, drift*lambda_c, group_cohesion, ...
                    y_init_coord(ag,rp), y_current_coord(ag,rp), cy, ...
                    dWy(ag,rp,t+burn_in_steps), board_size_y);

            elseif SIR_status(ag,rp,t) == 0
                % Susceptible
                if isempty(board{x,y,rp})
                    % No exposure here: base move + cohesion
                    SIR_status(ag,rp,t+1) = 0;
                    x_current_coord(ag,rp) = Next_Step_Cohesion( ...
                        diffusion, drift, group_cohesion, ...
                        x_init_coord(ag,rp), x_current_coord(ag,rp), cx, ...
                        dWx(ag,rp,t+burn_in_steps), board_size_x);

                    y_current_coord(ag,rp) = Next_Step_Cohesion( ...
                        diffusion, drift, group_cohesion, ...
                        y_init_coord(ag,rp), y_current_coord(ag,rp), cy, ...
                        dWy(ag,rp,t+burn_in_steps), board_size_y);

                else
                    % There is pathogen history: INCREMENTAL density + MULTIPLE draws (to match baseline)
                    pathogen_density = 0;
                    for ptg_time = board{x,y,rp}
                        pathogen_density = pathogen_density + ((1 - pathogen_decay_rate) ^ (t - ptg_time));

                        if rand() < 1 - exp(-infection_rate * pathogen_density)
                            % Becomes infected (same-step)
                            SIR_status(ag,rp,t+1) = 1;

                            % Move with infected-like params (λ-modified) + cohesion
                            x_current_coord(ag,rp) = Next_Step_Cohesion( ...
                                diffusion*lambda_d, drift*lambda_c, group_cohesion, ...
                                x_init_coord(ag,rp), x_current_coord(ag,rp), cx, ...
                                dWx(ag,rp,t+burn_in_steps), board_size_x);

                            y_current_coord(ag,rp) = Next_Step_Cohesion( ...
                                diffusion*lambda_d, drift*lambda_c, group_cohesion, ...
                                y_init_coord(ag,rp), y_current_coord(ag,rp), cy, ...
                                dWy(ag,rp,t+burn_in_steps), board_size_y);
                        else
                            % Stays susceptible (this iteration)
                            SIR_status(ag,rp,t+1) = 0;

                            % Move with susceptible params + cohesion
                            x_current_coord(ag,rp) = Next_Step_Cohesion( ...
                                diffusion, drift, group_cohesion, ...
                                x_init_coord(ag,rp), x_current_coord(ag,rp), cx, ...
                                dWx(ag,rp,t+burn_in_steps), board_size_x);

                            y_current_coord(ag,rp) = Next_Step_Cohesion( ...
                                diffusion, drift, group_cohesion, ...
                                y_init_coord(ag,rp), y_current_coord(ag,rp), cy, ...
                                dWy(ag,rp,t+burn_in_steps), board_size_y);
                        end
                    end

                end

            else
                % Recovered: base move + cohesion
                SIR_status(ag,rp,t+1) = 2;

                x_current_coord(ag,rp) = Next_Step_Cohesion( ...
                    diffusion, drift, group_cohesion, ...
                    x_init_coord(ag,rp), x_current_coord(ag,rp), cx, ...
                    dWx(ag,rp,t+burn_in_steps), board_size_x);

                y_current_coord(ag,rp) = Next_Step_Cohesion( ...
                    diffusion, drift, group_cohesion, ...
                    y_init_coord(ag,rp), y_current_coord(ag,rp), cy, ...
                    dWy(ag,rp,t+burn_in_steps), board_size_y);
            end
        end
    end
end

% -------- output folder (adds cohesion level) --------
out_sir = fullfile(output_folder, sprintf( ...
    'numAgents %d/numGroups %d/pathogenDecay %0.2f/recoveryRate %0.3f/infectionRate %0.2f/diffusionRate %0.2f/driftRate %0.2f/lambda_c %0.2f/lambda_d %0.2f/cohesion %0.3f', ...
    num_agents, num_groups, pathogen_decay_rate, recovery_rate, infection_rate, diffusion, drift, lambda_c, lambda_d, group_cohesion));

if ~exist(out_sir, 'dir')
    mkdir(out_sir);
end

csvwrite(fullfile(out_sir, 'MeanSirOUT.csv'),      txt_SIR);
csvwrite(fullfile(out_sir, 'StdSirOUT.csv'),       txt_std_SIR);
csvwrite(fullfile(out_sir, 'Peak.csv'),            txt_Peak);
csvwrite(fullfile(out_sir, 'Coordinates.csv'),     coord_record);
csvwrite(fullfile(out_sir, 'Raw_Last_Step.csv'),   txt_prevalence);

final_SIR = squeeze(SIR_status(:, :, end));
csvwrite(fullfile(out_sir, 'SIR_Status.csv'), final_SIR);
end


function data = Next_Sim_Step(diffusion, drift, init_coord, current_coord, dW, board_size, stp)
a = (init_coord - current_coord);
aa = (drift * a);
b = (diffusion * dW(:,:,stp));
data = current_coord + aa + b;
data(data < 1) = 1;
data(data >= board_size) = (board_size - 1);
end

% ---------- NEW: movement with cohesion term ----------
function next = Next_Step_Cohesion(diffusion, drift, cohesion, init_coord, current_coord, center_coord, dW, board_size)
% Hard-edge map (no torus)
next = current_coord ...
    + drift    * (init_coord   - current_coord) ...
    + diffusion* dW ...
    + cohesion * (center_coord - current_coord);

% Clamp to [1, board_size-1]
next(next < 1)           = 1;
next(next >= board_size) = (board_size - 1);
end
% -----------------------------------------------------


function data = Next_Step(diffusion, drift, init_coord, current_coord, dW, board_size)
data = current_coord + (drift * (init_coord - current_coord)) + (diffusion * dW);
data(data < 1) = 1;
data(data >= board_size) = (board_size - 1);
end

function data = create_random_3d_array(agents, repetitions, steps, burn_in_steps)
data = zeros(agents, repetitions, steps + burn_in_steps);
for i = 1:agents
    for j = 1:steps + burn_in_steps
        for k = 1:repetitions
            data(i, k, j) = randn();
        end
    end
end
end

function data = create_grid_world_array(board_size_x, board_size_y, num_repeats)
data = cell(board_size_x, board_size_y, num_repeats);
for i = 1:board_size_x
    for j = 1:board_size_y
        for k = 1:num_repeats
            data{i, j, k} = [];
        end
    end
end
end

function data = randomize_start(board_size, pop_size)
data = zeros(1, pop_size);
for i = 1:pop_size
    data(i) = rand() * (board_size - 1);
end
end

function data = randomize_group_start(board_size, num_agents, num_groups, num_repeats, num_per_group)
data = zeros(num_agents, num_repeats);
ct = 0;
aux = rand() * (board_size - 1);
for rep = 1:num_repeats
    for i = 1:num_agents
        if ct == num_per_group
            aux = rand() * (board_size - 1);
            ct  = 0;
        end
        ct = ct + 1;
        data(i, rep) = aux;
    end
end
end

function data = fixed_start(board_size, pop_size, initial_position)
data = ones(1, pop_size) * initial_position;
end

function data = infect_agents(num_infected_groups, num_infected_agents, num_per_group, num_agents, num_repeats)
data = zeros(num_agents, num_repeats);
for i = 1:num_infected_groups
    for j = 1:num_infected_agents
        data((i - 1) * num_per_group + j, :) = ones(1, num_repeats);
    end
end
end
