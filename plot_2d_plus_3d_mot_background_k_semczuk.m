title = '2d_plus_3d_mot_background_k_semczuk';
wd = run_example(title);
%% Load trajectories
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
z_center = 0.0015;  % z-coordinate of the centre

a = 0.0184;      % Semi-axis length along x
e = 0.94;
b = a*sqrt(1-e^2);      % Semi-axis length along y
c = b;      % Semi-axis length along z

% Define the grid and setup parameters
[x_grid, y_grid, z_grid] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.05, -0.05:0.001:0.05); % Create a 3D grid
% Set the parameters for the field calculation
direction = [1, 0, 0]; % Example direction
gradient_2D = 3.938; % Example gradient
detuning_2D = 32.99; % Detuning in MHz
gradient_3D = 8.5; % Example gradient
detuning_3D = detuning_2D; % Detuning in MHz
r_3d = 0.0075;
r_pipe = 0.00075;
r_push = 0.00195;
d = 0.35;
tolerance = 0.01; % Tolerance for matching B_L
u_B = 87.94; % MHz/T
h_bar = 1; % Reduced Planck's constant
r_L_2D = detuning_2D * h_bar / u_B * 1/gradient_2D; % Target magnetic field intensity
r_L_3D = detuning_3D * h_bar / u_B * 1/gradient_3D; % Target magnetic field intensity

% Create ellipsoid - Laser 2DMOT
[x, y, z] = ellipsoid(x_center, y_center, z_center, a, b, c, 50);
h = surf(x, y, z);
set(h, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Magnetic field 2DMOT
[x, y, z] = ellipsoid(x_center, y_center, z_center, 0.5 * r_L_2D, r_L_2D, r_L_2D, 50);
g = surf(x, y, z);
set(g, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Laser 3DMOT
[x, y, z] = ellipsoid(d, 0.0, 0.0, r_3d, r_3d, r_3d, 50);
k = surf(x, y, z);
set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Magnetic field 3DMOT
[x, y, z] = ellipsoid(d, 0.0, 0.0, r_L_3D, r_L_3D, 0.5 * r_L_3D, 50);
l = surf(x, y, z);
set(l, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Push beam
[x, y, z] = cylinder(r_pipe);
m = surf(0.45*z-0.05, y, x+z_center);
set(m, 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Small tube
[x, y, z] = cylinder(r_pipe);
n = surf((z)*0.02+0.05, y, x+z_center);
set(n, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Set axis limits and labels
% xlim([-0.05 0.40]);
% ylim([-0.02 0.02]);
% zlim([-0.02 0.02]);
fontsize(16,"points")
xlabel("X [m]");
ylabel("Y [m]");
zlabel("Z [m]");
grid on;
view([0 0]); 
axis equal;

hold off; % Stop adding plots to this axes

% Animation loop
i = 1;
while ishandle(f)
    % Update atom positions in the animation
    frame = position{i};
    set(atoms, 'XData', frame(:,1), 'YData', frame(:,2), 'ZData', frame(:,3));
    i = i + 1;
    
    xlim([ -0.05 0.4 ]);
    ylim([ -0.05 0.05 ]);
    zlim([ -0.05 0.05 ]);
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
ids_vel = [];
ids_pos = [];
ids_pos_5 = [];
ids_pos_7 = [];


for frame=output_vel'
    captured = sqrt(frame.vec(:,1).^2 + frame.vec(:,2).^2 + frame.vec(:,3).^2) < 0.3;
    ids_vel = unique([ids_vel; frame.id(captured)]);
end

for frame=output'
    captured_5 = (0.31 < frame.vec(:,1));
    ids_pos_5 = unique([ids_pos_5; frame.id(captured_5)]);
end

for frame=output'
    captured_7 = (0.39 < frame.vec(:,1));
    ids_pos_7 = unique([ids_pos_7; frame.id(captured_7)]);
end

ids_1 = intersect(ids_vel, ids_pos_5);
ids_2 = unique(setdiff(ids_1, ids_pos_7));

ids = ids_2;

% Note: the atoms will not have position vectors of equal length.
trajectories = cell(length(ids),1);

fprintf(mat2str(size(trajectories)), '\n');
fprintf(mat2str(size(ids_pos_7)), '\n');
fprintf(mat2str(size(ids_pos_5)), '\n');

for i=1:length(ids)
    id = ids(i);
    for frame=output'
        mask = frame.id == id;
        trajectories{i} = [ trajectories{i}; frame.vec(mask,:) ];
    end
end


%% Plot the trajectories
clf; set(gcf, 'Color', 'w');
for trajectory=trajectories'
    pos = trajectory{1};
    plot3(pos(:,1), pos(:,2), pos(:,3), '-k'); hold on;
end

axis equal;

% Close up on the source itself
xlim([ 0.30 0.40 ]);
ylim([ -0.05 0.05 ]);
zlim([ -0.05 0.05 ]);
xlabel("Position x [m]");
ylabel("Position y [m]");
zlabel("Position z [m]");

if isequal(title,'2d_plus_3d_mot_background_cs_semczuk')
    p = 2.0*10^(-9) * 133;
    R = 4.8;
    elip = 0.84;
    r = R*sqrt(1 - elip^2);
    S = (4 * 2*R * 2*r + 2 * 2*r * 2*r ) * 0.000001;
    T = 300.0;
    m = 133.0 * 1.660539 * 10^(-27);
    a = 1.0;
end

if isequal(title,'2d_plus_3d_mot_background_k_semczuk')
    p = 7.0*10^(-8) * 100;
    S = (4 * 36.8 * 12.56 + 2 * 12.56 * 12.56 ) * 0.000001;
    T = 316.15;
    m = 39.0 * 1.660539 * 10^(-27);
    a = 0.9326; %natural abundance
end

k = 1.380649 * 10^(-23);

trapping_rate = a*double(length(ids_pos_5))/double(10000000.0)*((S*(p))/(2.0*3.14*m*k*T)^(0.5));
c = length(ids_2)/double(21557.0);
display(trapping_rate);
display(double(length(ids_2))/double(10000000.0));

%%
clear plot;
a = 10^6 * 10^(-6);
b = 100;

x = (-300:0.5:300);
y_1 = a ./ (1 + (x * (1/b)).^2).^2;
y_2 = a * exp(-x.^2 ./ (2*(0.55*b)^2));

figure();
hold on;
plot(x, y_1,'Linewidth',2.0);
plot(x, y_2,'Linewidth',2.0);
hold off;
xlabel('{\it r} [mm]');
ylabel('{\it I} [m^{-2}s^{-1}]');
fontsize(16,'points');
grid on;
ax = gca; % Get current axes
ax.XTick = -300:60:300; % Set custom x-tick spacing
ax.YTick = 0:0.1:1; % Set custom y-tick spacing

% % Define the function to calculate the magnetic field
function B = calculate_field(pos, centre, gradient, direction)
    % Calculate the difference between position and centre
    delta = pos - centre;

    % Compute the component of delta along the direction
    z_comp = dot(delta, direction) * direction;

    % Compute the radial component (perpendicular to direction)
    r_comp = delta - z_comp;

    % Calculate the magnetic field
    B = gradient * (r_comp - 2 * z_comp); % Quadrupole formula
end