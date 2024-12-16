wd = run_example('3d_mot_oven_ag_uhlenberg');
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
f = figure(1);
clf; set(gcf, 'Color', 'w');
frame = position{1};
atoms = plot3(frame(:,1), frame(:,2), frame(:,3), '.k');
view([ 45 45 ]); axis equal;

xlim([ -0.005 0.005 ]);
ylim([ -0.5 0.005 ]);
zlim([ -0.005 0.005 ]);
xlabel("x");
ylabel("y");
zlabel("z");

% Define the ellipsoid parameters
x_center = 0;  % x-coordinate of the centre
y_center = 0;  % y-coordinate of the centre
z_center = 0;  % z-coordinate of the centre

R_2D_x = 0.034;      % Semi-axis length along x
R_2D_y = 0.008;      % Semi-axis length along y
R_2D_z = 0.008;      % Semi-axis length along z

R_3D = 0.0069;

r_push = 0.00075;

r_pipe = 0.00075;

% Define the grid and setup parameters
[x_grid, y_grid, z_grid] = meshgrid(-0.05:0.001:0.05, -0.05:0.001:0.05, -0.05:0.001:0.05); % Create a 3D grid
% Set the parameters for the field calculation
gradient_3D = 13.0; % Magnetic gradient in G/cm
detuning_3D = 9.918; % Detuning in MHz
tolerance = 0.01; % Tolerance for matching B_L
d = 0.5;
u_B = 1; % J/T
h_bar = 1; % Reduced Planck's constant J*s
r_Z_3D = detuning_3D * h_bar / u_B * 1/gradient_3D * 9.27 * 10^(-2); % Target magnetic field intensity

% Create ellipsoid - Laser 3DMOT
[x, y, z] = ellipsoid(d, 0.0, 0.0, R_3D, R_3D, R_3D, 50);
k = surf(x, y, z);
set(k, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create ellipsoid - Magnetic field 3DMOT
[x, y, z] = ellipsoid(d, 0.0, 0.0, r_Z_3D, r_Z_3D, 0.5 * r_Z_3D, 50);
l = surf(x, y, z);
set(l, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Push beam
[x, y, z] = cylinder(r_push);
m = surf(1*z-0.05, y, x+z_center);
set(m, 'FaceColor', [1 0.5 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

% Create cylinder - Small tube
[x, y, z] = cylinder(r_pipe);
n = surf((z)*0.02+0.05, y, x+z_center);
set(n, 'FaceColor', 'black', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Define the vertices of the rectangular box
x = ([0 1 1 0 0 1 1 0]-0.5)*0.1*5/5; % X-coordinates
y = ([0 0 1 1 0 0 1 1]-0.5)*0.1*2/5; % Y-coordinates
z = ([0 0 0 0 1 1 1 1]-0.5)*0.1*2/5+z_center; % Z-coordinates

% Define the edges of the box (pairs of vertices to connect)
edges = [
    1 2; 2 3; 3 4; 4 1; % Bottom face
    5 6; 6 7; 7 8; 8 5; % Top face
    1 5; 2 6; 3 7; 4 8  % Vertical edges
];

% Plot the edges of the box
for i = 1:size(edges, 1)
    plot3(x(edges(i, :)), y(edges(i, :)), z(edges(i, :)), 'g-', 'LineWidth', 2);
end

hold off;
axis equal;
grid on;

i=1;
while ishandle(f)
    frame = position{i};
    set(atoms, 'XData', frame(:,1), 'YData', frame(:,2), 'ZData', frame(:,3));
    i = i + 1;
    
    if i == length(position)
        i = 1;
        pause(0.1);
    end
    %title(sprintf('%d', i));
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
    captured = sqrt(frame.vec(:,1).^2 + frame.vec(:,2).^2 + frame.vec(:,3).^2) < 1.0;
    ids_vel = unique([ids_vel; frame.id(captured)]);
end

for frame=output'
    captured_5 = (-0.00338513 < frame.vec(:,2));
    ids_pos_5 = unique([ids_pos_5; frame.id(captured_5)]);
end

for frame=output'
    captured_7 = (0.0044 < frame.vec(:,2));
    ids_pos_7 = unique([ids_pos_7; frame.id(captured_7)]);
end

ids_2 = unique(setdiff(ids_pos_5, ids_pos_7));

ids = ids_vel;

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

if plot
    % Plot the trajectories
    clf; set(gcf, 'Color', 'w');
    for trajectory=trajectories'
        pos = trajectory{1};
        plot3(pos(:,1), pos(:,2), pos(:,3), '-k'); hold on;
    end
    
    axis equal;
    grid on;
    
    %% 
    % Close up on the source itself
    xlim([ -0.005 0.005 ]);
    ylim([ -0.005 0.005 ]);
    zlim([ -0.005 0.005 ]);
    xlabel("x");
    ylabel("y");
    zlabel("z");
end

p = 5.0*10^(-10) * 100;
S = (3.14*(0.564189)^2.0) * 0.000001;
T = 1200.0;
a = 0.482; %due to the limited abundance in the sample
m = 109.0 * 1.660539 * 10^(-27);
k = 1.380649 * 10^(-23);

A = 7.715;
B = -14935;
C = -0.2779;
D = -0.1701;
p_eq = 101325 * 10^(A + B/T + C*log10(T) + D*0.001*T);

evap_coff = 0.8;
%0.127 wynika z tego że przy losowaniu predkosci atomow bierzemy rozklad
%predkosci o temperaturze 300 a tak naprawde powinnismy 1200K
flux = 0.127 * 0.27 * a * double(length(ids))/double(2147483647.0)*((S*(p_eq - p))/(2.0*3.14*m*k*T)^(0.5));

display(flux);