% MAXIMUM SIMULATION TIME (days)
global tmax
tmax = 600;

% Initial conditions vector
x0 = zeros(2, 1);
x0(1) = 1.0; % Initial ROS concentration
x0(2) = 1.0; % Initial alpha-syn concentration

% Time span
tspan = 0:0.001:tmax;
opts = odeset('AbsTol', 1e-8);

% Simulate disease condition
disease_mode = true; 
[t_disease, x_disease] = ode23tb(@(t, x) f(t, x, disease_mode), tspan, x0, opts);

% Simulate non-disease condition
disease_mode = false;
[t_nondisease, x_nondisease] = ode23tb(@(t, x) f(t, x, disease_mode), tspan, x0, opts);

% Plot results
figure(1);

% Plot ROS concentration
subplot(2, 1, 1);
hold on
grid on
plot(t_disease, x_disease(:, 1), 'r', 'LineWidth', 2)         % disease conditions
plot(t_nondisease, x_nondisease(:, 1), 'b', 'LineWidth', 2);  % non-disease conditions
% vertical lines
%plot([300,300], [0, 100], 'Color', 'k', 'LineStyle','--') 
%plot([400,400], [0, 100], 'Color', 'k', 'LineStyle','--')
%plot([500,500], [0, 100], 'Color', 'k', 'LineStyle','--') 

legend('Disease (ROS)', 'Non-Disease (ROS)', 'Location','west');
title('ROS Concentration');
xlabel('Time (days)');
ylabel('ROS');

% Plot alpha-syn concentration
subplot(2, 1, 2);
hold on
grid on
plot(t_disease, x_disease(:, 2), 'r', 'LineWidth', 2);
plot(t_nondisease, x_nondisease(:, 2), 'b', 'LineWidth', 2);
%plot([500,500], [0, 100], 'Color', 'k', 'LineStyle','--')
%plot([300,300], [0, 100], 'Color', 'k', 'LineStyle','--') 
%plot([400,400], [0, 100], 'Color', 'k', 'LineStyle','--')
legend('Disease (alpha-syn)', 'Non-Disease (alpha-syn)', 'Location','west');
title('Alpha-Synuclein Concentration');
xlabel('Time (days)');
ylabel('Alpha-Syn');


% function to simulate the system of ODEs and set up the parameters
function xdot = f(t, x, disease_mode)
    global tmax
    % Compartment: id = Neuron, name = Neuron, constant
    compartment_Neuron = 1.0;    % VOLUME SCALING FACTOR
    % Parameters
    global_par_k1 = 17280.0;     % ROS GENERATION
    global_par_k2 = 17280.0;     % ROS REMOVAL
    global_par_k3 = 0.168;       % aSYN GENERATION
    global_par_k4 = 0.168;       % aSYN REMOVAL
    global_par_dalphasyn = 15.0; % increase in ROS release due to aSYN damage
    global_par_kalphasyn = 8.5;  % Hill kinetic constant for aSYN damage
    
    % set different S parameters for disease and normal states
    if disease_mode
        % global_par_S1 = piecewise(2.6, (t > 1) && (t < 300), 0); % simulate temporary input
        global_par_S1 = 2.6; % Oxidative stress (higher -> disease)
        % if t < 400
        %     global_par_S1 = 0 + 3*2.6*(t/tmax);
        % else
        %     global_par_S1 = 0;
        % end
        % global_par_S2 = piecewise(0.35, (t > 10) && (t < 300), 2); %
        global_par_S2 = 1; % Age-related anti-ox (higher -> no disease)
        % if t < 400
        %     global_par_S2 = 1 - 2*0.61*(t/tmax);
        % elseif t >= 400 && t < 500
        %     global_par_S2 = 1;
        % else
        %     global_par_S2 = 1.15;
        % end
        % global_par_S2 = 1 - 0.75*(t/tmax);
        global_par_S3 = 1; % Genetic predisposition (higher -> disease)
        global_par_S4 = 1; % Protein clearance (higher -> no disease)

    else % normal conditions
        global_par_S1 = 0; % Oxidative stress (higher -> disease)
        global_par_S2 = 1; % Age-related anti-ox (higher -> no disease)
        global_par_S3 = 1; % Genetic predisposition (higher -> disease)
        global_par_S4 = 1; % Protein clearance (higher -> no disease)
    end

    % REACTIONS
    reaction_ROS_1 = compartment_Neuron * V1(global_par_k1, global_par_S1, global_par_dalphasyn, x(2), global_par_kalphasyn);
    % SIMULATE REACTION EFFICIENCY INHIBITION
    % if t>100 && t<180
    %     reaction_ROS_1 = reaction_ROS_1 * 0.5; % inhibition magnitude
    % end

    % SIMULATE TWO DRUGS, DIFFERENT STRENGHTS
    % if t > 100 && t < 300
    %     if mod(t, 50) > 0 && mod(t, 50) < 15
    %         reaction_ROS_1 = reaction_ROS_1 * 0.5;
    %     end
    % end
    % if t > 500 && t < 700
    %     if mod(t, 50) > 0 && mod(t, 50) < 15
    %         reaction_ROS_1 = reaction_ROS_1 * 0.25;
    %     end
    % end
    reaction_ROS_2 = compartment_Neuron * V2(global_par_k2, x(1), global_par_S2);
    reaction_aSyn_1 = compartment_Neuron * V3(global_par_k3, x(1), global_par_S3);
    reaction_aSyn_2 = compartment_Neuron * V4(global_par_k4, x(2), global_par_S4);

    % ODEs
    xdot = zeros(2, 1);
    xdot(1) = (1 / compartment_Neuron) * ((1.0 * reaction_ROS_1) + (-1.0 * reaction_ROS_2));
    xdot(2) = (1 / compartment_Neuron) * ((1.0 * reaction_aSyn_1) + (-1.0 * reaction_aSyn_2));
end

% individual reactions
function z = V1(k1, Sx, d, S, k2)
    z = (k1 * (1 + Sx + d * (S / k2)^4 / (1 + (S / k2)^4)));
end

function z = V2(k, S, S2)
    z = (k * S * S2);
end

function z = V3(k, S, S3)
    z = (k * S * S3);
end

function z = V4(k, S, S4)
    z = (k * S * S4);
end

% set up time-dependant varying parameters
function z = piecewise(varargin)
    numArgs = nargin;
    result = 0;
    foundResult = 0;
    for k = 1:2:numArgs-1
        if varargin{k+1} == 1
            result = varargin{k};
            foundResult = 1;
            break;
        end
    end
    if foundResult == 0
        result = varargin{numArgs};
    end
    z = result;
end