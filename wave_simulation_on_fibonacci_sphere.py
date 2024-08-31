# Wave simulation on sphere surface (fibonacci sphere)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from mpl_toolkits.mplot3d import proj3d
from matplotlib.colors import Normalize
from scipy.spatial import KDTree


def fibonacci_sphere(num_points):
    global points
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
    for i in range(num_points):
        y = 1. - (i / float(num_points - 1)) * 2.  # y goes from 1 to -1
        radius = np.sqrt(1 - y ** 2)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))
    return np.array(points)


def update_arrow():
    global x_arw0, y_arw0, z_arw0, u_arw0, v_arw0, w_arw0
    global qvr_gaussian_center0
    global x_arw1, y_arw1, z_arw1, u_arw1, v_arw1, w_arw1
    global qvr_gaussian_center1
    # Arrow 0
    x_arw0 = np.sin(theta_gaussian0) * np.cos(phi_gaussian0)
    y_arw0 = np.sin(theta_gaussian0) * np.sin(phi_gaussian0)
    z_arw0 = np.cos(theta_gaussian0)
    u_arw0 = x_arw0 * 0.5
    v_arw0 = y_arw0 * 0.5
    w_arw0 = z_arw0 * 0.5
    qvr_gaussian_center0.remove()
    qvr_gaussian_center0 = ax0.quiver(x_arw0, y_arw0, z_arw0, u_arw0, v_arw0, w_arw0,
                                      length=1, color='red', normalize=False, label='Arrow0')
    # Arrow 1
    x_arw1 = np.sin(theta_gaussian1) * np.cos(phi_gaussian1)
    y_arw1 = np.sin(theta_gaussian1) * np.sin(phi_gaussian1)
    z_arw1 = np.cos(theta_gaussian1)
    u_arw1 = x_arw1 * 0.5
    v_arw1 = y_arw1 * 0.5
    w_arw1 = z_arw1 * 0.5
    qvr_gaussian_center1.remove()
    qvr_gaussian_center1 = ax0.quiver(x_arw1, y_arw1, z_arw1, u_arw1, v_arw1, w_arw1,
                                      length=1, color='blue', normalize=False, label='Arrow1')


def get_force(i_self):
    global tree
    dist, idx = tree.query(base_points[i_self], k=21)  # Set k=21 and find 21 points including itself.
    neighbors = idx[1:]     # 20 points excluding itself.
    dists = dist[1:]        # Distance of 20 points excluding itself.
    # force by magnitude
    force_by_magnitude = 0.
    for i in range(len(neighbors)):
        # effective_magnitude = (magnitude[neighbors[i]] - magnitude[i_self])
        effective_magnitude = (magnitude[neighbors[i]] - magnitude[i_self]) * (dists.mean() / dists[i]) ** 2.
        force_by_magnitude += k_spring * effective_magnitude
    return force_by_magnitude


def calc_wave():
    global points, xs, ys, zs
    global scat_sphere, magnitude, magnitude_buffer
    global velocity
    for i in range(number_of_points):
        force = get_force(i)
        a = force / mass
        velocity[i] = velocity[i] + a * 1.
        magnitude_buffer[i] = magnitude[i] + velocity[i] * 1.
    magnitude = magnitude_buffer.copy()
    points[:, 0] = (1. + magnitude) * base_points[:, 0]
    points[:, 1] = (1. + magnitude) * base_points[:, 1]
    points[:, 2] = (1. + magnitude) * base_points[:, 2]
    xs, ys, zs = zip(*points)


def update_scat():
    global scat_sphere
    scat_sphere.remove()
    scat_sphere = ax0.scatter(xs, ys, zs, c=magnitude, cmap=cmap_scat, s=size_scat, norm=norm)


def get_angle(point1, point2):
    magnitude1 = np.linalg.norm(point1)
    magnitude2 = np.linalg.norm(point2)
    dot_product = np.dot(point1, point2)
    cos_angle = dot_product / (magnitude1 * magnitude2)
    return np.arccos(cos_angle)


def apply_gaussian():
    global points, xs, ys, zs
    global scat_sphere, magnitude
    # print(index)
    # Gaussian 0 (red arrow)
    x0 = np.sin(theta_gaussian0) * np.cos(phi_gaussian0)
    y0 = np.sin(theta_gaussian0) * np.sin(phi_gaussian0)
    z0 = np.cos(theta_gaussian0)
    # Gaussian 1 (blue arrow)
    x1 = np.sin(theta_gaussian1) * np.cos(phi_gaussian1)
    y1 = np.sin(theta_gaussian1) * np.sin(phi_gaussian1)
    z1 = np.cos(theta_gaussian1)
    for i in range(number_of_points):
        # Gaussian 0 (red arrow)
        if var_chk0.get():
            angle0 = get_angle((base_points[i, 0], base_points[i, 1], base_points[i, 2]), (x0, y0, z0))
            gauss0 = scale_gaussian * (1. / (np.sqrt(2. * np.pi) * sigma) *
                                       np.exp(- (angle0 ** 2.) / (2. * sigma ** 2.)))
        else:
            gauss0 = 0.
        # Gaussian 1 (blue arrow)
        if var_chk1.get():
            angle1 = get_angle((base_points[i, 0], base_points[i, 1], base_points[i, 2]), (x1, y1, z1))
            gauss1 = scale_gaussian * (1. / (np.sqrt(2. * np.pi) * sigma) *
                                       np.exp(- (angle1 ** 2.) / (2. * sigma ** 2.)))
        else:
            gauss1 = 0.
        magnitude[i] = gauss0 + gauss1
        points[i, 0] = (1. + gauss0 + gauss1) * base_points[i, 0]
        points[i, 1] = (1. + gauss0 + gauss1) * base_points[i, 1]
        points[i, 2] = (1. + gauss0 + gauss1) * base_points[i, 2]
    xs, ys, zs = zip(*points)
    update_scat()


# Setter at tkinter
def set_theta0(value):
    global theta_gaussian_deg0, theta_gaussian0
    theta_gaussian_deg0 = int(value)
    theta_gaussian0 = np.deg2rad(theta_gaussian_deg0)
    update_arrow()


def set_phi0(value):
    global phi_gaussian_deg0, phi_gaussian0
    phi_gaussian_deg0 = int(value)
    phi_gaussian0 = np.deg2rad(phi_gaussian_deg0)
    update_arrow()


def set_theta1(value):
    global theta_gaussian_deg1, theta_gaussian1
    theta_gaussian_deg1 = int(value)
    theta_gaussian1 = np.deg2rad(theta_gaussian_deg1)
    update_arrow()


def set_phi1(value):
    global phi_gaussian_deg1, phi_gaussian1
    phi_gaussian_deg1 = int(value)
    phi_gaussian1 = np.deg2rad(phi_gaussian_deg1)
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


def set_k_spring(value):
    global k_spring
    reset()
    k_spring = float(value)


# Animation control
def step():
    global cnt
    global txt_step
    cnt += 1
    txt_step.set_text("Step=" + str(cnt))
    calc_wave()
    update_scat()


def reset():
    global is_play, cnt, txt_step
    global xs, ys, zs
    global velocity
    is_play = False
    cnt = 0
    txt_step.set_text("Step=" + str(cnt))
    velocity = np.zeros(number_of_points)
    for i in range(number_of_points):
        magnitude[i] = 0.
        points[i, 0] = base_points[i, 0]
        points[i, 1] = base_points[i, 1]
        points[i, 2] = base_points[i, 2]
    xs, ys, zs = zip(*points)
    update_scat()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt
    # global txt_step
    if is_play:
        step()
        # cnt += 1


# Global variables

# Animation control
cnt = 0
is_play = False

# Gaussian
sigma = 0.30
scale_gaussian = 0.05

theta_gaussian_deg0 = 0.
phi_gaussian_deg0 = 0.
theta_gaussian0 = np.deg2rad(theta_gaussian_deg0)
phi_gaussian0 = np.deg2rad(phi_gaussian_deg0)

theta_gaussian_deg1 = 90.
phi_gaussian_deg1 = - 90.
theta_gaussian1 = np.deg2rad(theta_gaussian_deg1)
phi_gaussian1 = np.deg2rad(phi_gaussian_deg1)

# Mass point and Spring constant
mass = 30.
k_spring = 2.
number_of_points = 4000

# Data structure
base_points = fibonacci_sphere(number_of_points)
points = base_points.copy()
xs, ys, zs = zip(*points)
magnitude = np.zeros(number_of_points)
magnitude_buffer = np.zeros(number_of_points)
velocity = np.zeros(number_of_points)

# Tree for search points
tree = KDTree(base_points)

# Generate figure and axes
title_tk = 'Wave on sphere surface (fibonacci sphere)'
title_ax0 = title_tk

range_xyz = 1.1
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
# ax0.set_facecolor('lightblue')
ax0.xaxis.pane.set_facecolor('darkgray')
ax0.yaxis.pane.set_facecolor('darkgray')
ax0.zaxis.pane.set_facecolor('darkgray')

# Generate items
# Text items
txt_step = ax0.text2D(x_min0, y_max0, "Step=" + str(cnt))
xz, yz, _ = proj3d.proj_transform(x_min0, y_max0, z_max0, ax0.get_proj())
txt_step.set_position((xz, yz))

# Plot items
size_scat = 8
cmap_scat = 'seismic'
norm = Normalize(vmin=-0.05, vmax=0.05)
scat_sphere = ax0.scatter(xs, ys, zs, c=magnitude, cmap=cmap_scat, s=size_scat, norm=norm)

# Arrow
x_arw0 = np.sin(theta_gaussian0) * np.cos(phi_gaussian0)
y_arw0 = np.sin(theta_gaussian0) * np.sin(phi_gaussian0)
z_arw0 = np.cos(theta_gaussian0)
u_arw0 = x_arw0 * 0.5
v_arw0 = y_arw0 * 0.5
w_arw0 = z_arw0 * 0.5
qvr_gaussian_center0 = ax0.quiver(x_arw0, y_arw0, z_arw0, u_arw0, v_arw0, w_arw0,
                                  length=1, color='red', normalize=False, label='Arrow0')

x_arw1 = np.sin(theta_gaussian1) * np.cos(phi_gaussian1)
y_arw1 = np.sin(theta_gaussian1) * np.sin(phi_gaussian1)
z_arw1 = np.cos(theta_gaussian1)
u_arw1 = x_arw1 * 0.5
v_arw1 = y_arw1 * 0.5
w_arw1 = z_arw1 * 0.5
qvr_gaussian_center1 = ax0.quiver(x_arw1, y_arw1, z_arw1, u_arw1, v_arw1, w_arw1,
                                  length=1, color='blue', normalize=False, label='Arrow1')

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
lbl_ks = tk.Label(frm_parameters, text="k(spring constant)")
lbl_ks.pack(side='left')
var_ks = tk.StringVar(root)  # variable for spinbox-value
var_ks.set(str(k_spring))  # Initial value
spn_ks = tk.Spinbox(
    frm_parameters, textvariable=var_ks, format="%.2f", from_=1., to=10.0, increment=1.,
    command=lambda: set_k_spring(float(var_ks.get())), width=5
)
spn_ks.pack(side='left')
lbl_mass = tk.Label(frm_parameters, text="Point mass")
lbl_mass.pack(side='left')
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(str(mass))  # Initial value
spn_mass = tk.Spinbox(
    frm_parameters, textvariable=var_mass, format="%.2f", from_=1., to=50., increment=1.,
    command=lambda: set_mass(float(var_mass.get())), width=5
)
spn_mass.pack(side='left')

# gaussian curve center (arrow position)
# Arrow 0
frm_arw0 = ttk.Labelframe(root, relief="ridge", text="Gaussian center0(red arrow)", labelanchor="n", width=100)
frm_arw0.pack(side='left')

var_chk0 = tk.BooleanVar(frm_arw0)
var_chk0.set(True)
chk0 = tk.Checkbutton(frm_arw0, text="", variable=var_chk0)
chk0.pack(side='left')

lbl_theta0 = tk.Label(frm_arw0, text='Theta:')
lbl_theta0.pack(side='left')
var_theta0 = tk.StringVar(root)
var_theta0.set(str(theta_gaussian_deg0))
spn_theta0 = tk.Spinbox(
    frm_arw0, textvariable=var_theta0, format='%.0f', from_=-180, to=180, increment=1,
    command=lambda: set_theta0(var_theta0.get()), width=8
    )
spn_theta0.pack(side='left')

lbl_phi0 = tk.Label(frm_arw0, text='Phi:')
lbl_phi0.pack(side='left')
var_phi0 = tk.StringVar(root)
var_phi0.set(str(phi_gaussian_deg0))
spn_phi0 = tk.Spinbox(
    frm_arw0, textvariable=var_phi0, format='%.0f', from_=-360, to=360, increment=1,
    command=lambda: set_phi0(var_phi0.get()), width=8
    )
spn_phi0.pack(side='left')

# Arrow 1
frm_arw1 = ttk.Labelframe(root, relief="ridge", text="Gaussian center1(blue arrow)", labelanchor="n", width=100)
frm_arw1.pack(side='left')

var_chk1 = tk.BooleanVar(frm_arw1)
var_chk1.set(False)
chk1 = tk.Checkbutton(frm_arw1, text="", variable=var_chk1)
chk1.pack(side='left')

lbl_theta1 = tk.Label(frm_arw1, text='Theta:')
lbl_theta1.pack(side='left')
var_theta1 = tk.StringVar(root)
var_theta1.set(str(theta_gaussian_deg1))
spn_theta1 = tk.Spinbox(
    frm_arw1, textvariable=var_theta1, format='%.0f', from_=-180, to=180, increment=1,
    command=lambda: set_theta1(var_theta1.get()), width=8
    )
spn_theta1.pack(side='left')

lbl_phi1 = tk.Label(frm_arw1, text='Phi:')
lbl_phi1.pack(side='left')
var_phi1 = tk.StringVar(root)
var_phi1.set(str(phi_gaussian_deg1))
spn_phi1 = tk.Spinbox(
    frm_arw1, textvariable=var_phi1, format='%.0f', from_=-360, to=360, increment=1,
    command=lambda: set_phi1(var_phi1.get()), width=8
    )
spn_phi1.pack(side='left')

# gaussian curve
frm_gaussian = ttk.Labelframe(root, relief="ridge", text="Gaussian", labelanchor="n", width=100)
frm_gaussian.pack(side='left')

lbl_sigma = tk.Label(frm_gaussian, text="Sigma")
lbl_sigma.pack(side='left')
var_sigma = tk.StringVar(root)  # variable for spinbox-value
var_sigma.set(str(sigma))  # Initial value
spn_sigma = tk.Spinbox(
    frm_gaussian, textvariable=var_sigma, format="%.2f", from_=0.01, to=1.0, increment=0.01,
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
