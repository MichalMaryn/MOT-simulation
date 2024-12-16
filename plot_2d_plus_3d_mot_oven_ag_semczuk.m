%wd = run_example('3d_mot_oven_ag_semczuk_effusion_cell');
wd = run_example('2d_plus_3d_mot_oven_ag_semczuk');
%% Load trajectories and plot results.
%
plot = true;

output = read_output(fullfile(wd, 'pos.txt'));
output_vel = read_output(fullfile(wd, 'vel.txt'));
position = {output.vec};
velocity = {output_vel.vec};
fprintf('Positions loaded.\n')

p_txt = fullfile('atomecs','pos.txt');
v_txt = fullfile('atomecs','vel.txt');
p_csv = fullfile('atomecs','pos.csv');
v_csv = fullfile('atomecs','vel.csv');

YourPos = readtable(p_txt, Delimiter=[",", ': ', ')', '(']);
YourVel = readtable(v_txt, Delimiter=[",", ': ', ')', '(']);

A = find(sum(ismissing(YourPos),2)==5);

YourPos(:,1) = [];
YourPos(:,2) = [];

YourVel(:,1) = [];
YourVel(:,2) = [];

YourPos = renamevars(YourPos,["Var2","Var4","Var5","Var6"],["Key1","PosX","PosY", "PosZ"]);
YourVel = renamevars(YourVel,["Var2","Var4","Var5","Var6"],["Key1","VelX","VelY", "VelZ"]);

writetable(YourPos, p_csv);
writetable(YourVel, v_csv);
%%

% Prepare the figure and plot atoms
f = figure(1);
clf; set(gcf, 'Color', 'w');
frame = position{1};
atoms = plot3(frame(:,1), frame(:,2), frame(:,3), '.k'); % Plot initial frame
hold on;  % Allow adding additional plots to the same axes

% Define the ellipsoid parameters
x_center = 0;  % x-coordinate of the centre
y_center = 0;  % y-coordinate of the centre
z_center = 0.002;  % z-coordinate of the centre

a = 0.01244;      % Semi-axis length along x
e = 0.56;
b = a*sqrt(1-e^2);      % Semi-axis length along y
c = b;      % Semi-axis length along z

% Define the grid and setup parameters
[x_grid, y_grid, z_grid] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.05, -0.05:0.001:0.05); % Create a 3D grid
% Set the parameters for the field calculation
direction = [1, 0, 0]; % Example direction
gradient_2D = 12.885; % Example gradient
detuning_2D = 63.72; % Detuning in MHz
gradient_3D = 8.5; % Example gradient
detuning_3D = detuning_2D; % Detuning in MHz
r_3d = 0.0075;
r_pipe = 0.00075;
r_push = 0.00162;
d = 0.21;
tolerance = 0.01; % Tolerance for matching B_L
u_B = 87.94; % MHz/T
h_bar = 1; % Reduced Planck's constant
r_L_2D = detuning_2D * h_bar / u_B * 1/gradient_2D; % Target magnetic field intensity
r_L_3D = detuning_3D * h_bar / u_B * 1/gradient_3D; % Target magnetic field intensity

% Create ellipsoid - Laser 2DMOT
[x, y, z] = ellipsoid(x_center, y_center, 0, a, b, c, 50);
h = surf(x, y, z);
set(h, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Magnetic field 2DMOT
[x, y, z] = ellipsoid(x_center, y_center, 0, 0.5 * r_L_2D, r_L_2D, r_L_2D, 50);
g = surf(x, y, z);
set(g, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Laser 3DMOT
[x, y, z] = ellipsoid(d, 0.0, -z_center, r_3d, r_3d, r_3d, 50);
k = surf(x, y, z);
set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Magnetic field 3DMOT
[x, y, z] = ellipsoid(d, 0.0, -z_center, r_L_3D, r_L_3D, 0.5 * r_L_3D, 50);
l = surf(x, y, z);
set(l, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Push beam
[x, y, z] = cylinder(r_pipe);
m = surf(0.45*z-0.05, y, x);
set(m, 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Small tube
[x, y, z] = cylinder(r_pipe);
n = surf((z)*0.01+0.02, y, x);
set(n, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Set axis limits and labels
xlim([-0.05 0.40]);
ylim([-0.02 0.02]);
zlim([-0.02 0.02]);
xlabel("Position x [m]");
ylabel("Position y [m]");
zlabel("Position z [m]");
grid on;
view([0 0]); 
% axis equal;

hold off; % Stop adding plots to this axes

% Animation loop
i = 1;
while ishandle(f)
    % Update atom positions in the animation
    frame = position{i};
    set(atoms, 'XData', frame(:,1), 'YData', frame(:,2), 'ZData', frame(:,3));
    i = i + 1;
    
    % Reset to the first frame when reaching the end
    if i == length(position)
        i = 1;
        pause(0.2);
    end
    pause(0.1);
end
%%
% Identify the atoms that are pushed out of the 2D MOT, and plot their
% trajectories.

ids = [];
ids_2d = [];
ids_pos_5 = [];

for frame=output_vel'
    captured = sqrt(frame.vec(:,1).^2 + frame.vec(:,2).^2 + frame.vec(:,3).^2) < 0.2;
    ids = unique([ids; frame.id(captured)]);
end

for frame=output'
    captured_2d = (sqrt(frame.vec(:,1).^2 + (1.0/(1-0.56^2))*frame.vec(:,2).^2 + (1.0/(1.0-0.56^2))*frame.vec(:,3).^2) < 0.012885.^2.0);
    ids_2d = unique([ids_2d; frame.id(captured_2d)]);
end

for frame=output'
    captured_5 = (0.2 < frame.vec(:,1));
    ids_pos_5 = unique([ids_pos_5; frame.id(captured_5)]);
end

% Note: the atoms will not have position vectors of equal length.
trajectories = cell(length(ids),1);

fprintf(mat2str(size(trajectories)), '\n');
fprintf(mat2str(size(ids_2d)), '\n');

for i=1:length(ids)
    id = ids(i);
    for frame=output'
        mask = frame.id == id;
        trajectories{i} = [ trajectories{i}; frame.vec(mask,:) ];
    end
end

if plot
    % Plot the trajectories
    clf; set(gcf, 'Color', 'w');
    for trajectory=trajectories'
        pos = trajectory{1};
        plot3(pos(:,1), pos(:,2), pos(:,3), '-k'); hold on;
    end
    
    axis equal;
    
    %% 
    % Close up on the source itself
    xlim([ -0.03 0.2175 ]);
    ylim([ -0.016 0.016 ]);
    zlim([ -0.020 0.016 ]);
    xlabel("x");
    ylabel("y");
    zlabel("z");
end

% with the increase of temperature we have less atoms of speed below capture velocity
% and the flux decreases based on the equation but when the temperature
% decreases the evaporation preassure decreases significantly
T = 1173.15; 
A = 7.715;
B = -14935;
C = -0.2779;
D = -0.1701;
p_eq = 101325 * 10^(A + B/T + C*log10(T) + D*0.001*T);

p = 5.0*10^(-10) * 100;
S = 1.0 * 18.278 * 0.000001;
m = 109.0 * 1.660539 * 10^(-27);
k = 1.380649 * 10^(-23);
eff_v = 0.9; % this is 0.9 when we evaporate something from a solid blok of metal, this is not the case since it is both an effusion cell
F_N = 0.333; %this is a mistake take the 0.9 it is true for 850K
a = 0.4865 ;
flux = a*F_N * double(length(ids))/double(10000000.0)*((S*(p_eq - p))/(2.0*3.14*m*k*T)^(0.5));

display(F_N * double(length(ids))/double(10000000.0)*((S*(p_eq - p))/(2.0*3.14*m*k*T)^(0.5)));
display(F_N * double(length(ids_pos_5))/double(10000000.0)*((S*(p_eq - p))/(2.0*3.14*m*k*T)^(0.5)));

display(flux);