import numpy as np
from numpy import diff
from datetime import datetime
import csv
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
import ast

atom = 'Cesium'
if atom == 'Potassium':
    radius_push = 0.75 * 0.001
    power_push = 0.1 * 0.001
    detuning_push = -4.978 * 6.035

    radius_2d = 18.38 * 0.001
    power_2d = 48.803 * 0.001
    detuning_2d = -4.978 * 6.035
    ellipticity = 0.943

    radius_3d = 3.0 * 0.001
    power_3d = 35.729 * 0.001
    detuning_3d = -35.003

    quadrupole_gradient_2d = 5.018
    quadrupole_gradient_3d = 18.397

    distance_oven_2d = 19.9 * 0.001
    distance_2d_3d = 350.0 * 0.001
    shift_x = 19.9 * 0.001
    shift_z = 1.0 * 0.001
    thick = 1.0 * 0.001

    temperature = 316.15

    nb0 = 5000000.0
    velocity_cap = 70.0

    Grav = 'G'
    t = 6000

if atom == 'Potassium_catani':
    detuning_push = -5.2 * 6.035;
    power_push = 0.006;
    radius_push = 0.75e-3;
    wavelenght = 767.0e-9;

    detuning_2d = -5.8 * 6.035;
    power_2d = 0.08;
    radius_2d = 14.1e-3;
    ellipticity = 0.943;

    radius_3d = 0.0
    power_3d = 0.0
    detuning_3d = 0.0

    quadrupole_gradient_2d = 17.0
    quadrupole_gradient_3d = 0.0

    distance_oven_2d = 0.0 * 0.001
    distance_2d_3d = 0.0 * 0.001
    shift_x = 0.0 * 0.001
    shift_z = 1.0 * 0.001
    thick = 1.0 * 0.001

    temperature = 316.15

    nb0 = 10000000.0
    velocity_cap = 70.0

    Grav = 'G'
    t = 600

if atom == 'Cesium':
    radius_push = 1.40285259522301 * 0.001
    power_push = 0.885678766695696 * 0.001
    detuning_push = -5.03 * 5.22

    radius_2d = 4.78294334304728 * 0.001
    power_2d = 98.9067599877347 * 0.001
    detuning_2d = -5.03 * 5.22
    ellipticity = 0.842044819918941

    radius_3d = 7.5 * 0.001
    power_3d = 35.493 * 0.001
    detuning_3d = -5.03 * 5.22

    quadrupole_gradient_2d = 8.386
    quadrupole_gradient_3d = 6.369

    distance_oven_2d = 0.0 * 0.001
    distance_2d_3d = 350.0 * 0.001
    shift_x = 0.0 * 0.001
    shift_z = 1.49608569809914 * 0.001
    thick = 1.0 * 0.001

    temperature = 300.0

    nb0 = 5000000.0
    velocity_cap = 50.0

    Grav = 'G'
    t = 6000

if atom == 'Cesium_lam':
    radius_push = 1.5e-3
    power_push = 3.0e-3
    detuning_push = -1.9 * 5.22

    radius_2d = 34.0e-3
    power_2d = 150.0e-3
    detuning_2d = -1.9 * 5.22
    ellipticity = 0.943

    radius_3d = 6.9e-3
    power_3d = 150.0e-3
    detuning_3d = -1.9 * 5.22

    quadrupole_gradient_2d = 13.0
    quadrupole_gradient_3d = 8.9

    distance_oven_2d = 0.0 * 0.001
    distance_2d_3d = 610.0 * 0.001
    shift_x = 0.0 * 0.001
    shift_z = 3.0 * 0.001
    thick = 1.0 * 0.001

    temperature = 300.0

    nb0 = 40000000.0
    velocity_cap = 70.0

    Grav = 'G'
    t = 6000

if atom == 'Silver':
    radius_push = 1.62 * 0.001
    power_push = 0.09 * 0.001
    detuning_push = -2.72 * 23.4

    radius_2d = 12.44 * 0.001
    power_2d = 249.0 * 0.001
    detuning_2d = -2.72 * 23.4
    ellipticity = 0.56

    radius_3d = 7.5 * 0.001
    power_3d = 50.0 * 0.001
    detuning_3d = -2.72 * 23.4

    quadrupole_gradient_2d = 12.278
    quadrupole_gradient_3d = 5.0

    distance_oven_2d = 15.0 * 0.001
    distance_2d_3d = 210.0 * 0.001
    shift_x = 0.0 * 0.001
    shift_z = 0.0 * 0.001
    thick = 0.0 * 0.001

    temperature = 900.0 + 273.15

    nb0 = 10000000.0
    velocity_cap = 90.0

    Grav = 0
    t = 8000

small_radius = (radius_2d**2.0 - ellipticity**2.0 * radius_2d**2.0)**0.5

# def read_param_2D():

#     path = 'atomecs/pos.csv'
#     path_v = 'atomecs/vel.csv'

#     df = pd.read_csv(path, dtype='float', delimiter=',')
#     df_v = pd.read_csv(path_v, dtype='float', delimiter=',')

#     df.drop(df.loc[df['Key1']==0.0].index, inplace=True)
#     df_v.drop(df_v.loc[df_v['Key1']==0.0].index, inplace=True)

#     # df_id = df.isna().sum(axis=1)
#     # steps = np.where(np.array(df_id.to_list())==3)[0]

#     # ids = set(df.loc[((df_v.iloc[:,1]<60.0) & (df.iloc[:,1]>0.0146)),'Key1'].to_list())
#     ids = set(df.loc[((df.loc[:,'PosX']>0.03)),'Key1'].to_list())
#     ids_all = df.isna().loc[:, 'PosX'].index[df.isna().loc[:, 'PosX']].to_numpy()
#     nb1 = len(ids)

#     print(nb1)

#     posis = []
#     velos = []

#     cooling_t = []

#     last_vel = []

#     for i in ids:

#         pos = df.loc[df.iloc[:,0] == i].drop('Key1', axis=1).to_numpy()
#         cool_t = np.where(pos[:,0]>0.03)[0][0]
#         cooling_t.append(cool_t)
#         pos = np.pad(pos, ((0,t-len(pos)),(0,0)), 'constant', constant_values=np.nan)

#         vel = df_v.loc[df_v.iloc[:,0] == i].drop('Key1', axis=1).to_numpy()
#         last_vel.append(np.sqrt(vel[-1][0]**2 + vel[-1][1]**2 + vel[-1][2]**2))
#         vel = np.pad(vel, ((0,t-len(vel)),(0,0)), 'constant', constant_values=np.nan)

#         if np.exp(-60*np.array(cool_t)/100000) > 0.5: 
#             posis.append(pos)
#             velos.append(vel)

#     velos = np.array(velos)
#     posis = np.array(posis)

#     return posis, velos, cooling_t, [nb0, nb1]

def read_param_3D():

    path = 'atomecs/pos.csv'
    path_v = 'atomecs/vel.csv'

    df = pd.read_csv(path, dtype='float', delimiter=',')
    df_v = pd.read_csv(path_v, dtype='float', delimiter=',')

    df.drop(df.loc[df['Key1']==0.0].index, inplace=True)
    df_v.drop(df_v.loc[df_v['Key1']==0.0].index, inplace=True)

    # ids_all = df.isna().loc[:, 'PosX'].index[df.isna().loc[:, 'PosX']].to_numpy()

    # ids = df.loc[(ids_all[len(ids_all)-1]+1):,'Key1']

    ids = np.unique(df.loc[df['PosX']>0.05, 'Key1'].to_numpy())

    nb1 = len(np.unique(df.loc[df['PosX']>0.02, 'Key1'].to_numpy()))

    time_index = np.unique(df.loc[df['PosX']>0.02, 'Key1'].index)

    posis = []
    initial_velos = []
    velos_small = []
    vel_x_at_pos = []
    set_pos = 0.1
    capture_velocity = 0

    for i in ids:
        # t_2dmot = 0
        pos = df.loc[df['Key1']==i].drop(['Key1'], axis = 1).to_numpy()
        # posis.append(pos)
        vel = df_v.loc[df_v['Key1']==i].drop(['Key1'], axis = 1).to_numpy()
        final_vel = np.sqrt(vel[-1][0]**2 + vel[-1][1]**2 + vel[-1][2]**2)
        if final_vel < 1.0 and pos[-1][0] > 0.56 and pos[-1][0] < 0.66 and atom == 'Cesium_lam':
            posis.append(pos)
            velos_small.append(np.rot90(vel, k=1, axes=(0, 1)))
            # vel_Y_Z = np.array(np.sqrt(vel[:, 1]**2 + vel[:, 2]**2))
            # t_2dmot = (np.argmax(vel_Y_Z<0.5))
            if np.argmax(pos[:,0]>set_pos) != 0:
                # if ((np.argmax(pos[:,0]>set_pos))-t_2dmot) <= 0:
                #     t_2dmot = 0
                vel_x_at_pos.append(set_pos / ((np.argmax(pos[:,0]>set_pos))-0)*10000) # Extract VelX
        if final_vel < 1.0 and pos[-1][0] > 0.31 and pos[-1][0] < 0.39 and (atom == 'Cesium' or atom == 'Potassium'):
            posis.append(pos)
            velos_small.append(np.rot90(vel, k=1, axes=(0, 1)))
            vel_Y_Z = np.array(np.sqrt(vel[:, 1]**2 + vel[:, 2]**2))
            t_2dmot = (np.argmax(vel_Y_Z<0.5))
            if np.argmax(pos[:,0]>set_pos) != 0:
                if ((np.argmax(pos[:,0]>set_pos))-t_2dmot) <= 0:
                    t_2dmot = 0
                vel_x_at_pos.append(set_pos / ((np.argmax(pos[:,0]>set_pos))-t_2dmot)*10000) # Extract VelX
        if pos[-1][0] > 0.02 and atom == 'Potassium_catani':
            velos_small.append(np.rot90(vel, k=1, axes=(0, 1)))
            velos_small.append(np.rot90(vel, k=1, axes=(0, 1)))
            # vel_Y_Z = np.array(np.sqrt(vel[:, 1]**2 + vel[:, 2]**2))
            # t_2dmot = (np.argmax(vel_Y_Z<0.5))
            if np.argmax(pos[:,0]>set_pos) != 0:
                # if ((np.argmax(pos[:,0]>set_pos))-t_2dmot) <= 0:
                #     t_2dmot = 0
                vel_x_at_pos.append(set_pos / ((np.argmax(pos[:,0]>set_pos))-0)*10000) # Extract VelX
        if final_vel < 1.0 and pos[-1][0] > 0.2026 and pos[-1][0] < 0.2174 and atom == 'Silver':
            velos_small.append(np.rot90(vel, k=1, axes=(0, 1)))
            vel_Y_Z = np.array(np.sqrt(vel[:, 1]**2 + vel[:, 2]**2))
            t_2dmot = (np.argmax(vel_Y_Z<0.5))
            if np.argmax(pos[:,0]>set_pos) != 0:
                if ((np.argmax(pos[:,0]>set_pos))-t_2dmot) <= 0:
                    t_2dmot = 0
                vel_x_at_pos.append(set_pos / ((np.argmax(pos[:,0]>set_pos))-t_2dmot)*10000) 
        initial_velos.append(np.sqrt(vel[0][0]**2 + vel[0][1]**2 + vel[0][2]**2))
        # condition = np.isclose(pos[:, 0], set_pos, atol=0.01)  # Adjust tolerance if needed

    print(vel_x_at_pos)
    vel_x_at_pos = np.array(vel_x_at_pos)
    # vel_x_at_pos = vel_x_at_pos[vel_x_at_pos > 3]
    # vel_x_at_pos = vel_x_at_pos[vel_x_at_pos < 20]

    # Plot histogram
    # plt.figure()
    # counts, bins, _ = plt.hist(vel_x_at_pos, bins=10, edgecolor='black', alpha=0.6, density=True, label="Histogram")
    # plt.close()

    # Fit Gaussian
    # mu, sigma = norm.fit(vel_x_at_pos)

    # Generate Gaussian curve
    # x = np.linspace(bins[0], bins[-1], 1000)
    # pdf = norm.pdf(x, mu, sigma)

    # Define capture criteria
    # percentile = 95  # Define the desired capture percentage
    # capture_velocity = norm.ppf(percentile / 100, loc=mu, scale=sigma)

    # Plot Gaussian fit
    # plt.plot(x, pdf, 'r-', label=f"Gaussian Fit\n$\\mu$ = {mu:.2f}, $\\sigma$ = {sigma:.2f}")
    # plt.title(f"Histogram and Gaussian Fit of X-Velocities at Position {set_pos}")
    # plt.xlabel("X-Velocity")
    # plt.ylabel("Probability Density")
    # plt.legend()
    # plt.show()

    #print(nb1)
    eff_2D = nb1/nb0
    eff_3D = len(velos_small)/nb1

    #print(len(velos_small))

    # plt.figure()
    # plt.hist(initial_velos, bins=10)
    # plt.show()

    return np.array(posis) , velos_small, [eff_2D, eff_3D], time_index, capture_velocity

# def load_t_2D():

#     path = 'atomecs/pos.txt'
#     path_v = 'atomecs/vel.txt'

#     df = pd.read_csv(path, dtype='string', delimiter=' ')
#     df_v = pd.read_csv(path_v, dtype='string', delimiter=' ')

#     all_ids = df.head(30000)[['step-100,']].to_numpy()

#     trapped = []

#     for i in all_ids:
#         i = i[0]
#         if i.split('-')[0] != 'step':
#             pos = df[df.iloc[:,0] == i].iloc[:,1].to_numpy()
#             pos = np.array([ast.literal_eval(i) for i in pos])
#             if any(pos[:,0] > 0.5):
#                 print(i)
#                 vel = df_v[df_v.iloc[:,0] == i].iloc[:,1].to_numpy()
#                 vel = np.array([ast.literal_eval(i) for i in vel])
#                 vel = np.sqrt(vel[:,1]**2 + vel[:,2]**2)[-400:]
#                 when_not_0 = np.where(((vel - np.roll(vel, 1)) < 0))[0]
#                 when_not_0 = when_not_0[len(when_not_0)-1]
#                 if when_not_0 < 199: trapped.append(when_not_0)

#     return np.mean(trapped)/10

# def plot3D(ids, pos):

#     fig = plt.figure()
#     ax = plt.axes(projection='3d')

#     for i in range(len(ids)):
#         ax.plot(pos[i,:,0], pos[i,:,1],zs=pos[i,:,2])
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     plt.show()

# def velocity_distr(posis, velos, posX):

#     vel_profile = []

#     for j in range(len(posis)):
#         for i in range(int(t/10)):
#             if posis[j, i, 0] > posX:
#                 vel_profile.append(np.sqrt(abs(velos[j, i, 0]**2 + velos[j, i, 1]**2 + velos[j, i, 2]**2)))

#     plt.figure()
#     plt.hist(vel_profile, bins=10)
#     plt.show()

#     return vel_profile

# def save_to_file(nb):

#     start = datetime.now()

#     with open('result_log.csv', 'a') as f_object:
#         writer = csv.writer(f_object)
#         # writer.writerow(['Date', "Atom", 'NbAtoms', '2DRad [m]', '2DEli', '2DPower [W]', '2DDet', '2DZPos [m]', 'PPower [W]', 'PRad [m]', 'PDet', "2D-3D dist [m]", '3DRad [m]', '3DEli', '3DPower [W]', '3DDet [m]', '3DZPos [m]', 'Grav [m/s]', 'Time [ms]', 'Eff [e-5]'])
#         writer.writerow([start, atom, nb[0], radius_push, power_push, detuning_push, radius_2d, power_2d, detuning_2d, ellipticity, radius_3d, power_3d, detuning_3d, distance_oven_2d, distance_2d_3d,shift_x, shift_z, thick, Grav, t, (nb[1]/nb[0])*100000])
#         f_object.close()

cooling_t = []

posis, velos, eff, time_index, capture_velocity = read_param_3D()
# posis, velos, time, nb = read_param_2D()

print(np.sum(time_index)/(100000*len(time_index)))

# a = velocity_distr(posis, velos, 0.1)

if atom == 'Potassium':
    atom_u = 39.0
    p = 2.1*10**(-7) * 100 #preasure in N/m^2 
    S = 2 * (small_radius*2.0)**2 + 4 * (radius_2d*2.0)*(small_radius*2.0) # sourface in meters squared

if atom == 'Potassium_catani':
    atom_u = 39.0
    p = 7.0*10**(-8) * 100 #preasure in N/m^2 
    S = 2 * (small_radius*2.0)**2 + 4 * (radius_2d*2.0)*(small_radius*2.0) # sourface in meters squared

if atom == 'Cesium' or atom == 'Cesium_lam':
    atom_u = 133.0
    p = 1*10**(-9) * 133 #preasure in N/m^2 
    S = 2 * (small_radius*2.0)**2 + 4 * (radius_2d*2.0)*(small_radius*2.0) # sourface in meters squared

if atom == 'Silver':
    atom_u = 109.0
    p = 5.0*10**(-10) * 100
    # S = 8.0 * ( 0.25 * np.pi * 3.175 ** 2.0 + np.pi * 3.175 * ( 3.175 - 0.127 ) ) * 10 ** -6.0
    S = 1.0 * 18.2 * 0.000001

m = atom_u * 1.660539 * 10**(-27) #in kg/mol * mol = kg
k = 1.380649 * 10**(-23)
flux = eff[0]*((S*p)/np.sqrt(2.0*np.pi*m*k*temperature))

print('Atom: ', atom)
print('Atomic flux from a 2D MOT: ', flux/10000000000, ' * 10^(10) atom/s')
print('Efficiency of trapping for a 3D MOT: ', eff[1])
print('Trapping rate in 3D MOT: ', eff[1]*flux/10000000, ' * 10^(7) atom/s')
print('Capture velocity of a 2D MOT: ', capture_velocity, ' m/s')

# save_to_file(nb)

# print(nb[1]/nb[0]*100000)
        
# print(load_t_2D())

# print(velos[0][0][0])

# velxyz = np.array(np.sqrt(velos[0,0]**2 + velos[0,1]**2 + velos[0,2]**2))
# velyz = np.array(np.sqrt(velos[:,:,1]**2 + velos[:,:,2]**2))

# fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3)

# for i in velxyz:
#     ax2.plot(np.arange(len(i)),i)

# for i in velyz:
#     ax3.plot(np.arange(len(i)),i)

# for i in velos[:,:,0]:
#     ax4.plot(np.arange(len(i)),i)

# for i in velos[:,:,1]:
#     ax5.plot(np.arange(len(i)),i)

# for i in velos[:,:,2]:
#     ax6.plot(np.arange(len(i)),i)

# ax2.set_title('XYZ vel')
# ax3.set_title('YZ vel')
# ax4.set_title('X vel')
# ax5.set_title('Y vel')
# ax6.set_title('Z vel')
# plt.show()

abs_v = []
plt.figure()
for i in np.arange(len(velos)):
    if i%20 == 0:
        v = np.sqrt(velos[i][0]**2 + velos[i][1]**2 + velos[i][2]**2)
        abs_v.append(v)
        plt.plot(np.arange(0,len(velos[i][0])/10.0,0.1),v)
        # plt.plot(np.resize(posis[i], (3,600))[0],v)
plt.xlabel('Time [ms]')
plt.ylabel('Velocity [m/s]')
plt.title('Absolute Velocity of Trapped Atoms')
plt.xlim(0, 45.0)
plt.savefig('velo_time_cs_semczuk.png', dpi=1000)
plt.show()

df = pd.DataFrame(abs_v)
df.to_csv('velocities.csv', mode='a', index=False, header=False)
# what to do: correct the number of all atoms by the velocity cap, velocity distr and loading rate with a one group of atoms by measruing how many atoms are in 2D MOT check the definition of a loading rate