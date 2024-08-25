# Wave simulation on sphere surface test
# Note; Due to the conversion from Cartesian coordinates to polar coordinates,
# the mesh near the poles and the equator is uneven,
# causing the wave simulation to not work well.

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def update_arrow():
    global x_arw, y_arw, z_arw, u_arw,v_arw,w_arw
    global qvr_gaussian_center
    x_arw = np.sin(theta_gaussian) * np.cos(phi_gaussian)
    y_arw = np.sin(theta_gaussian) * np.sin(phi_gaussian)
    z_arw = np.cos(theta_gaussian)
    u_arw = x_arw * 0.5
    v_arw = y_arw * 0.5
    w_arw = z_arw * 0.5
    qvr_gaussian_center.remove()
    qvr_gaussian_center = ax0.quiver(x_arw, y_arw, z_arw, u_arw, v_arw, w_arw,
                                     length=1, color='red', normalize=False, label='Arrow')


def get_opposite(num):
    return int((num + len(phi_linspace) / 2) % len(phi_linspace))


def get_force(i, j):
    # Phi direction
    if i == 0:
        dz_i = - (magnitude[i + 1][j] - magnitude[i][j]) - (magnitude[-1][j] - magnitude[i][j])
    elif i == len(phi_linspace) - 1:
        dz_i = - (magnitude[i - 1][j] - magnitude[i][j]) - (magnitude[0][j] - magnitude[i][j])
    else:
        dz_i = - (magnitude[i - 1][j] - magnitude[i][j]) - (magnitude[i + 1][j] - magnitude[i][j])
    force_i = - k * dz_i
    # Theta direction
    i_opposite = get_opposite(i)
    if j == 0:
        dz_j = - (magnitude[i][j + 1] - magnitude[i][j]) # - (magnitude[i_opposite][j] - magnitude[i][j])
    elif j == len(theta_linspace) - 1:
        dz_j = - (magnitude[i][j - 1] - magnitude[i][j]) # - (magnitude[i_opposite][j] - magnitude[i][j])
    else:
        dz_j = - (magnitude[i][j - 1] - magnitude[i][j]) - (magnitude[i][j + 1] - magnitude[i][j])
    force_j = - k * dz_j
    force = force_i + force_j
    return force


def update_sphere():
    # global wireframe
    global plt_sphere
    global x, y, z
    x = (1. + magnitude) * np.sin(theta) * np.cos(phi)
    y = (1. + magnitude) * np.sin(theta) * np.sin(phi)
    z = (1. + magnitude) * np.cos(theta)
    plt_sphere.remove()
    plt_sphere = ax0.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=plt.cm.seismic(norm(magnitude)),
                                  edgecolor='black', linewidth=0.1, alpha=1)


def apply_gaussian():
    global magnitude
    '''
    magnitude = phi * 0. + theta * 0.
    magnitude = 1 / (2. * np.pi * sigma ** 2.) * np.exp(-((phi - phi_gaussian) ** 2. + (theta - theta_gaussian) ** 2.) / (2. * sigma ** 2.))
    magnitude *= scale_gaussian
    '''
    magnitude_center = 1. / (2. * np.pi * sigma ** 2.) * np.exp(-((phi - np.pi) ** 2. + (theta - theta_gaussian) ** 2.)
                                                                / (2. * sigma ** 2.))
    index_roll = int(((len(phi_linspace) - 0) / (2 * np.pi)) * (phi_gaussian - np.pi))
    magnitude = np.roll(magnitude_center, index_roll, axis=0)
    magnitude *= scale_gaussian
    update_sphere()


# Setter at tkinter
def set_theta(value):
    global theta_gaussian_deg, theta_gaussian
    theta_gaussian_deg = int(value)
    theta_gaussian = np.deg2rad(theta_gaussian_deg)
    update_arrow()


def set_phi(value):
    global phi_gaussian_deg, phi_gaussian
    phi_gaussian_deg = int(value)
    phi_gaussian = np.deg2rad(phi_gaussian_deg)
    update_arrow()


def set_sigma(value):
    global sigma
    sigma = float(value)


def set_scale_gaussian(value):
    global scale_gaussian
    scale_gaussian = float(value)


def set_mass(value):
    global mass
    reset()
    mass = float(value)


def set_k(value):
    global k
    reset()
    k = float(value)


# Animation control
def step():
    global cnt
    cnt += 1


def reset():
    global is_play, cnt, magnitude, velocity
    is_play = False
    cnt = 0
    txt_step.set_text("Step=" + str(cnt))
    magnitude = phi * 0. + theta * 0.
    velocity = magnitude * 0.
    update_sphere()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt, txt_step, magnitude, magnitude_buffer, velocity
    if is_play:
        txt_step.set_text("Step=" + str(cnt))
        for i in range(len(phi_linspace)):
            for j in range(len(theta_linspace)):
                force = get_force(i, j)
                a = force / mass
                velocity[i][j] = velocity[i][j] + a * 1.
                magnitude_buffer[i][j] = magnitude[i][j] + velocity[i][j] * 1.
        magnitude = magnitude_buffer.copy()
        velocity_row_means = np.mean(velocity, axis=1)
        magnitude_row_means = np.mean(magnitude, axis=1)
        for ix in range(len(phi_linspace)):
            velocity[ix][0] = velocity_row_means[0]
            velocity[ix][-1] = velocity_row_means[-1]
            magnitude[ix][0] = magnitude_row_means[0]
            magnitude[ix][-1] = magnitude_row_means[-1]
        update_sphere()
        cnt += 1


# Global variables

# Animation control
cnt = 0
is_play = False

# Gaussian
sigma = 0.1
scale_gaussian = 0.05
theta_gaussian_deg = 90.
phi_gaussian_deg = - 90.
theta_gaussian = np.deg2rad(theta_gaussian_deg)
phi_gaussian = np.deg2rad(phi_gaussian_deg)

# Mass point and Spring constant
mass = 20.
k = 10.

# Data structure
theta_linspace = np.linspace(0.01, np.pi - 0.01, 100)
# theta_linspace = np.arccos(1. - 2. * np.linspace(0., 1., 100))
phi_linspace = np.linspace(0., 2. * np.pi, 200)
theta, phi = np.meshgrid(theta_linspace, phi_linspace)
magnitude = theta * 0. + phi * 0.

x = (1. + magnitude) * np.sin(theta) * np.cos(phi)
y = (1. + magnitude) * np.sin(theta) * np.sin(phi)
z = (1. + magnitude) * np.cos(theta)

magnitude_buffer = magnitude * 0.
velocity = magnitude * 0.

# Generate figure and axes
title_tk = 'Wave on sphere surface'
title_ax0 = title_tk

range_xyz = 1.2
x_min0 = - range_xyz
x_max0 = range_xyz
y_min0 = - range_xyz
y_max0 = range_xyz
z_min0 = - range_xyz
z_max0 = range_xyz

fig = Figure()
ax0 = fig.add_subplot(111, projection='3d')
ax0.set_box_aspect((1, 1, 1))

ax0.set_xlim(x_min0, x_max0)
ax0.set_ylim(y_min0, y_max0)
ax0.set_zlim(z_min0, z_max0)
ax0.set_title(title_ax0)
ax0.set_xlabel('x')
ax0.set_ylabel('y')
ax0.set_zlabel('z')
ax0.grid()

# Generate items
# Text items
txt_step = ax0.text2D(x_min0, y_max0, "Step=" + str(cnt))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))

# Plot items
# wireframe = ax0.plot_wireframe(x, y, z, rstride=4, cstride=4)
norm = Normalize(vmin=-0.2, vmax=0.2)
plt_sphere = ax0.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=plt.cm.seismic(norm(magnitude)),
                              edgecolor='black', linewidth=0.1, alpha=1)
# Arrow
x_arw = np.sin(theta_gaussian) * np.cos(phi_gaussian)
y_arw = np.sin(theta_gaussian) * np.sin(phi_gaussian)
z_arw = np.cos(theta_gaussian)
u_arw = x_arw * 0.5
v_arw = y_arw * 0.5
w_arw = z_arw * 0.5
qvr_gaussian_center = ax0.quiver(x_arw, y_arw, z_arw, u_arw, v_arw, w_arw,
                                 length=1, color='red', normalize=False, label='Arrow')

# Tkinter
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Animation
frm_anim = ttk.Labelframe(root, relief='ridge', text='Animation', labelanchor='n')
frm_anim.pack(side='left', fill=tk.Y)
btn_play = tk.Button(frm_anim, text='Play/Pause', command=switch)
btn_play.pack(side='left')
btn_step = tk.Button(frm_anim, text='Step', command=step)
btn_step.pack(side='left')
btn_reset = tk.Button(frm_anim, text='Reset', command=reset)
btn_reset.pack(side='left')

# Parameters
frm_parameters = ttk.Labelframe(root, relief="ridge", text="Parameters", labelanchor="n", width=100)
frm_parameters.pack(side='left')
lbl_k = tk.Label(frm_parameters, text="k(spring constant)")
lbl_k.pack(side='left')
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(str(k))  # Initial value
spn_k = tk.Spinbox(
    frm_parameters, textvariable=var_k, format="%.2f", from_=1., to=10.0, increment=1.,
    command=lambda: set_k(float(var_k.get())), width=5
)
spn_k.pack(side='left')
lbl_mass = tk.Label(frm_parameters, text="Point mass")
lbl_mass.pack(side='left')
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(str(mass))  # Initial value
spn_mass = tk.Spinbox(
    frm_parameters, textvariable=var_mass, format="%.2f", from_=20., to=50., increment=1.,
    command=lambda: set_mass(float(var_mass.get())), width=5
)
spn_mass.pack(side='left')

# gaussian curve
frm_gaussian = ttk.Labelframe(root, relief="ridge", text="Gaussian", labelanchor="n", width=100)
frm_gaussian.pack(side='left')


lbl_theta = tk.Label(frm_gaussian, text='Theta:')
lbl_theta.pack(side='left')
var_theta = tk.StringVar(root)
var_theta.set(str(theta_gaussian_deg))
spn_theta = tk.Spinbox(
    frm_gaussian, textvariable=var_theta, format='%.0f', from_=-180, to=180, increment=1,
    command=lambda: set_theta(var_theta.get()), width=8
    )
spn_theta.pack(side='left')

lbl_phi = tk.Label(frm_gaussian, text='Phi:')
lbl_phi.pack(side='left')
var_phi = tk.StringVar(root)
var_phi.set(str(phi_gaussian_deg))
spn_phi = tk.Spinbox(
    frm_gaussian, textvariable=var_phi, format='%.0f', from_=-360, to=360, increment=1,
    command=lambda: set_phi(var_phi.get()), width=8
    )
spn_phi.pack(side='left')

lbl_sigma = tk.Label(frm_gaussian, text="Sigma")
lbl_sigma.pack(side='left')
var_sigma = tk.StringVar(root)  # variable for spinbox-value
var_sigma.set(str(sigma))  # Initial value
spn_sigma = tk.Spinbox(
    frm_gaussian, textvariable=var_sigma, format="%.2f", from_=0.1, to=1.0, increment=0.1,
    command=lambda: set_sigma(float(var_sigma.get())), width=5
)
spn_sigma.pack(side='left')

lbl_scale = tk.Label(frm_gaussian, text="Scale:")
lbl_scale.pack(side='left')
var_scale = tk.StringVar(root)  # variable for spinbox-value
var_scale.set(str(scale_gaussian))  # Initial value
spn_scale = tk.Spinbox(
    frm_gaussian, textvariable=var_scale, format="%.3f", from_=0.001, to=2.0, increment=0.001,
    command=lambda: set_scale_gaussian(var_scale.get()), width=5
)
spn_scale.pack(side='left')

btn_gaussian = tk.Button(frm_gaussian, text="Apply", command=lambda: apply_gaussian())
btn_gaussian.pack(side='left')

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=50, save_count=100)
root.mainloop()
