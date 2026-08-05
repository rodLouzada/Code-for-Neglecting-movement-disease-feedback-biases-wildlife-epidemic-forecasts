function Simulation_Framework()

    rng(42, 'multFibonacci');

    sir_mode = 1;  % Change this variable to '0' (zero) to make an SI simulation
    
    output_folder = 'path\to\folder\Torus MultiLambda Z Score SIR Status';% 
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
    
    fprintf('Starting Simulation\n');
    
    simulation(num_agents, ...
               num_groups, ...
               board_size_x, ...
               board_size_y, ...
               num_infected_agents, ...
               num_infected_groups, ...
               recovery_rate, ...
               pathogen_decay, ...
               pathogen_deposition, ...
               num_steps, ...
               diffusion, ...
               drift, ...
               time_step, ...
               num_repeats, ...
               infection_rate, ...
               infected_movement_impact_c,infected_movement_impact_d, ...
               sir_mode, ...
               output_folder);
                                        
                                        
end

function simulation(num_agents, num_groups, board_size_x, board_size_y, num_infected_agents, num_infected_groups, ...
               recovery_rate, pathogen_decay_rate, pathogen_deposition, num_steps, diffusion, drift, time_step, ...
               num_repeats, infection_rate, lambda_c, lambda_d, sir_mode, output_folder)
    
    num_per_group = num_agents / num_groups;
    burn_in_steps = 10;
    
    txt_SIR = zeros(num_steps, 3);
    txt_Peak = zeros(num_repeats, 2);
    txt_std_SIR = zeros(num_steps, 3);
    txt_prevalence = zeros(num_repeats, 3);
    
    SIR_status = zeros(num_agents, num_repeats, num_steps+1);
    
    coord_record = zeros(num_agents* num_steps,5); %step,ID, X, Y, SIR
    
    dWx = create_random_3d_array(num_agents, num_repeats, num_steps, burn_in_steps);
    dWy = create_random_3d_array(num_agents, num_repeats, num_steps, burn_in_steps);
    
    board = create_grid_world_array(board_size_x, board_size_y, num_repeats);
                  
    x_init_coord = randomize_group_start(board_size_x, num_agents, num_groups, num_repeats, num_per_group);    
    y_init_coord = randomize_group_start(board_size_y, num_agents, num_groups, num_repeats, num_per_group);
    
    x_current_coord = x_init_coord;
    y_current_coord = y_init_coord;
    
    % Burn-in step
    for bis = 1:burn_in_steps
        x_current_coord = Next_Sim_Step(diffusion, drift, x_init_coord, x_current_coord, dWx, board_size_x, bis);
        y_current_coord = Next_Sim_Step(diffusion, drift, y_init_coord, y_current_coord, dWy, board_size_y, bis);
    end
    
    % Infect the initial agents 
    SIR_status(:,:,1) = infect_agents(num_infected_groups, num_infected_agents, num_per_group, num_agents, num_repeats);
    
    % For every timestep
    for t = 1:num_steps
        txt_SIR(t,1) = sum(sum(SIR_status(:,:,t) == 0)) / num_repeats;
        txt_SIR(t,2) = sum(sum(SIR_status(:,:,t) == 1)) / num_repeats;
        txt_SIR(t,3) = sum(sum(SIR_status(:,:,t) == 2)) / num_repeats;

         % Compute the standard deviation
        txt_std_SIR(t,1) = std(sum(SIR_status(:,:,t) == 0, 1));
        txt_std_SIR(t,2) = std(sum(SIR_status(:,:,t) == 1, 1));
        txt_std_SIR(t,3) = std(sum(SIR_status(:,:,t) == 2, 1));

        % Store prevalence at the last timestep
        if t == num_steps
            txt_prevalence(:,1) = sum(SIR_status(:,:,t) == 0, 1); % S at last step per repetition
            txt_prevalence(:,2) = sum(SIR_status(:,:,t) == 1, 1); % I at last step per repetition
            txt_prevalence(:,3) = sum(SIR_status(:,:,t) == 2, 1); % R at last step per repetition
        end

        
        for rp = 1:num_repeats
            
            if sum(SIR_status(:,rp,t) == 1) > txt_Peak(rp,1)
                txt_Peak(rp,1) = sum(SIR_status(:,rp,t));
                txt_Peak(rp,2) = t;
            end
            
            for ag = 1:num_agents
                x = round(x_current_coord(ag,rp));
                y = round(y_current_coord(ag,rp));
                
                if(rp==1)
                        coord_record((ag-1)*num_steps + (t),:) = [t,ag,x_current_coord(ag,rp),y_current_coord(ag,rp),SIR_status(ag,rp,t)]; %step,ID, X, Y, SIR
                end

                if SIR_status(ag,rp,t) == 1
                    
                    % Deposit pathogen
                    board{x,y,rp} = [board{x,y,rp}, t];
                    
                    % Try to recover
                    if rand() < recovery_rate
                        SIR_status(ag,rp,t+1) = 2;
                    else
                        SIR_status(ag,rp,t+1) = 1;
                    end
                    
                    % Move to the next spot
                    x_current_coord(ag,rp) = Next_Step(diffusion*lambda_d, drift*lambda_c, x_init_coord(ag,rp), x_current_coord(ag,rp), dWx(ag,rp,t+burn_in_steps), board_size_x);
                    y_current_coord(ag,rp) = Next_Step(diffusion*lambda_d, drift*lambda_c, y_init_coord(ag,rp), y_current_coord(ag,rp), dWy(ag,rp,t+burn_in_steps), board_size_y);
                    
                elseif SIR_status(ag,rp,t) == 0
                    
                    % If not in an infected place
                    if isempty(board{x,y,rp})
                        x_current_coord(ag,rp) = Next_Step(diffusion, drift, x_init_coord(ag,rp), x_current_coord(ag,rp), dWx(ag,rp,t+burn_in_steps), board_size_x);
                        y_current_coord(ag,rp) = Next_Step(diffusion, drift, y_init_coord(ag,rp), y_current_coord(ag,rp), dWy(ag,rp,t+burn_in_steps), board_size_y);
                        SIR_status(ag,rp,t+1) = 0;
                    else
                        pathogen_density = 0;
                        
                        for ptg_time = board{x,y,rp}
                            pathogen_density = pathogen_density + ((1 - pathogen_decay_rate) ^ ( t- ptg_time)); %(t- ptg_time)
                            
                            if rand() < 1 - exp(-infection_rate * pathogen_density)
                                SIR_status(ag,rp,t+1) = 1;
                                x_current_coord(ag,rp) = Next_Step(diffusion*lambda_d, drift*lambda_c, x_init_coord(ag,rp), x_current_coord(ag,rp), dWx(ag,rp,t+burn_in_steps), board_size_x);
                                y_current_coord(ag,rp) = Next_Step(diffusion*lambda_d, drift*lambda_c, y_init_coord(ag,rp), y_current_coord(ag,rp), dWy(ag,rp,t+burn_in_steps), board_size_y);
                            else
                                SIR_status(ag,rp,t+1) = 0;
                                x_current_coord(ag,rp) = Next_Step(diffusion, drift, x_init_coord(ag,rp), x_current_coord(ag,rp), dWx(ag,rp,t+burn_in_steps), board_size_x);
                                y_current_coord(ag,rp) = Next_Step(diffusion, drift, y_init_coord(ag,rp), y_current_coord(ag,rp), dWy(ag,rp,t+burn_in_steps), board_size_y);
                            end
                        end
                    end
                    
                elseif SIR_status(ag,rp,t) == 2
                    SIR_status(ag,rp,t+1) = 2;
                    
                    x_current_coord(ag,rp) = Next_Step(diffusion, drift, x_init_coord(ag,rp), x_current_coord(ag,rp), dWx(ag,rp,t+burn_in_steps), board_size_x);
                    y_current_coord(ag,rp) = Next_Step(diffusion, drift, y_init_coord(ag,rp), y_current_coord(ag,rp), dWy(ag,rp,t+burn_in_steps), board_size_y);
                end
            end
        end
    end
    
    out_sir = fullfile(output_folder, sprintf('numAgents %d/numGroups %d/pathogenDecay %0.2f/recoveryRate %0.3f/infectionRate %0.2f/diffusionRate %0.2f/driftRate %0.2f/lambda_c %0.2f/lambda_d %0.2f', ...
        num_agents, num_groups, pathogen_decay_rate, recovery_rate, infection_rate, diffusion, drift, lambda_c,lambda_d));
    if ~exist(out_sir, 'dir')
        mkdir(out_sir);
    end
    
    csvwrite(fullfile(out_sir, 'MeanSirOUT.csv'), txt_SIR);
    csvwrite(fullfile(out_sir, 'StdSirOUT.csv'), txt_std_SIR);
    csvwrite(fullfile(out_sir, 'Peak.csv'), txt_Peak);
    csvwrite(fullfile(out_sir, 'Coordinates.csv'), coord_record);
    csvwrite(fullfile(out_sir, 'Raw_Last_Step.csv'), txt_prevalence);
    final_SIR = squeeze( SIR_status(:, :, end) );
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


function data = Next_Step(diffusion, drift, init_coord, current_coord, dW, board_size)
    % Calculate the difference in coordinates considering periodic boundaries
    a = init_coord - current_coord;
    % Adjust the drift component to account for wraparound distance
    % The mod operation is adjusted to handle periodic boundary conditions
    wrap_a = mod(a + board_size/2, board_size) - board_size/2;
    % Compute drift and diffusion components
    aa = drift * wrap_a;
    b = diffusion * dW;
    % Update the coordinates
    data = current_coord + aa + b;
    % Apply wraparound logic to ensure coordinates stay within bounds
    % The modulo operation is adjusted to ensure proper periodic boundary handling
    data = mod(data - 1, board_size) + 1;
    % Ensure the coordinates are integers within bounds
    data = round(data);
    % Handle potential out-of-bounds errors due to rounding
    data(data < 1) = 1;
    data(data > board_size) = board_size;
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
                ct = 0;
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
