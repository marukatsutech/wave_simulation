# Wave simulation 3D
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)


def get_force(i, j):
    force = 0.
    if var_boundary_cond.get() == 1:
        if i == 0:
            dz_i = 0.
        elif i == 1:
            dz_i = - (0. - Z[i][j]) - (Z[i + 1][j] - Z[i][j])
        elif i == len(x) - 2:
            dz_i = - (Z[i - 1][j] - Z[i][j]) - (0. - Z[i][j])
        elif i == len(x) - 1:
            dz_i = 0.
        else:
            dz_i = - (Z[i - 1][j] - Z[i][j]) - (Z[i + 1][j] - Z[i][j])
        force_i = - k * dz_i
        if j == 0:
            dz_j = 0.
        elif j == 1:
            dz_j = - (0. - Z[i][j]) - (Z[i][j + 1] - Z[i][j])
        elif j == len(y) - 2:
            dz_j = - (Z[i][j - 1] - Z[i][j]) - (0. - Z[i][j])
        elif j == len(y) - 1:
            dz_j = 0.
        else:
            dz_j = - (Z[i][j - 1] - Z[i][j]) - (Z[i][j + 1] - Z[i][j])
        force_j = - k * dz_j
    else:
        if i == 0:
            dz_i = - (Z[i + 1][j] - Z[i][j])
        elif i == len(x) - 1:
            dz_i = - (Z[i - 1][j] - Z[i][j])
        else:
            dz_i = - (Z[i - 1][j] - Z[i][j]) - (Z[i + 1][j] - Z[i][j])
        force_i = - k * dz_i
        if j == 0:
            dz_j = - (Z[i][j + 1] - Z[i][j])
        elif j == len(y) - 1:
            dz_j = - (Z[i][j - 1] - Z[i][j])
        else:
            dz_j = - (Z[i][j - 1] - Z[i][j]) - (Z[i][j + 1] - Z[i][j])
        force_j = - k * dz_j
    force = force_i + force_j
    return force


def redraw_wireframe():
    global wireframe
    wireframe.remove()
    wireframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)


def impact_gaussian():
    global Z
    Z = X * 0. + Y * 0.
    x_start = distance_gaussian / 2.
    if num_gaussian == 1:
        x_step = 0.
    else:
        x_step = distance_gaussian / (num_gaussian - 1)
    for i in range(num_gaussian):
        Z = Z + 1 / (2. * np.pi * sigma ** 2.) * np.exp(-((X - x_start + i * x_step) ** 2. + Y ** 2.) / (2. * sigma ** 2.))
    redraw_wireframe()


def change_distance_gaussian(value):
    global distance_gaussian
    distance_gaussian = value
    impact_gaussian()


def change_num_gaussian(value):
    global num_gaussian
    num_gaussian = value
    impact_gaussian()


def change_sigma(value):
    global sigma
    sigma = value
    impact_gaussian()


def change_mass(value):
    global mass
    reset()
    mass = value


def change_k(value):
    global k
    reset()
    k = value


def set_axis():
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    ax.set_title('Wave simulation 3D')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_box_aspect((1, 1, 0.5))
    ax.grid()


def reset():
    global is_play, cnt, Z, wireframe, V
    is_play = False
    cnt = 0
    txt_step.set_text("Step=" + str(cnt))
    Z = X * 0. + Y * 0.
    V = Z * 0.
    redraw_wireframe()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt, txt_step, Z, wireframe, Z_buffer, V
    if is_play:
        txt_step.set_text("Step=" + str(cnt))
        for i in range(len(x)):
            for j in range(len(y)):
                force = get_force(i, j)
                a = force / mass
                V[i][j] = V[i][j] + a * 1.
                Z_buffer[i][j] = Z[i][j] + V[i][j] * 1.
        Z = Z_buffer.copy()
        redraw_wireframe()
        cnt += 1


# Global variables
x_min = -4
x_max = 4.
y_min = -4.
y_max = 4.
z_min = -4.
z_max = 4.

# Parameters
sigma = 0.2
mass = 20.
k = 10.

# Mesh grid
x = np.arange(x_min, x_max, 0.1)
y = np.arange(y_min, y_max, 0.1)
X, Y = np.meshgrid(x, y)
Z = X * 0. + Y * 0.

Z_buffer = Z
V = Z * 0.

# For UI
is_play = False
cnt = 0

num_gaussian = 1
distance_gaussian = 0.

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect((2, 2, 1))
set_axis()

# Generate items
txt_step = ax.text(x_min, y_max, z_max, "Step=" + str(0))
wireframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)

# Tkinter
root = tk.Tk()
root.title("Wave simulation 3D")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Parameters
frm_parameters = ttk.Labelframe(root, relief="ridge", text="Parameters", labelanchor="n", width=100)
frm_parameters.pack(side='left')
lbl_k = tk.Label(frm_parameters, text="k(spring constant)")
lbl_k.pack(side='left')
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(str(k))  # Initial value
spn_k = tk.Spinbox(
    frm_parameters, textvariable=var_k, format="%.2f", from_=1., to=10.0, increment=1.,
    command=lambda: change_k(float(var_k.get())), width=5
)
spn_k.pack(side='left')
lbl_mass = tk.Label(frm_parameters, text="Point mass")
lbl_mass.pack(side='left')
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(str(mass))  # Initial value
spn_mass = tk.Spinbox(
    frm_parameters, textvariable=var_mass, format="%.2f", from_=20., to=50., increment=1.,
    command=lambda: change_mass(float(var_mass.get())), width=5
)
spn_mass.pack(side='left')
# Boundary condition
frm_boundary_cond = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
frm_boundary_cond.pack(side='left')
var_boundary_cond = tk.IntVar(root)
rdb_fxe = tk.Radiobutton(frm_boundary_cond, text="Fixed end, ", value=1, var=var_boundary_cond)
rdb_fxe.pack(side='left')
rdb_fre = tk.Radiobutton(frm_boundary_cond, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack(side='left')
var_boundary_cond.set(1)

# gaussian curve
frm_gaussian = ttk.Labelframe(root, relief="ridge", text="Gaussian", labelanchor="n", width=100)
frm_gaussian.pack(side='left')
lbl_sigma = tk.Label(frm_gaussian, text="Sigma")
lbl_sigma.pack(side='left')
var_sigma = tk.StringVar(root)  # variable for spinbox-value
var_sigma.set(str(sigma))  # Initial value
spn_sigma = tk.Spinbox(
    frm_gaussian, textvariable=var_sigma, format="%.2f", from_=0.1, to=1.0, increment=0.1,
    command=lambda: change_sigma(float(var_sigma.get())), width=5
)
spn_sigma.pack(side='left')
lbl_num_gauss = tk.Label(frm_gaussian, text="Number")
lbl_num_gauss.pack(side='left')
var_num_gauss = tk.StringVar(root)  # variable for spinbox-value
var_num_gauss.set(str(num_gaussian))  # Initial value
spn_num_gauss = tk.Spinbox(
    frm_gaussian, textvariable=var_num_gauss, from_=1, to=6, increment=1,
    command=lambda: change_num_gaussian(int(var_num_gauss.get())), width=5
)
spn_num_gauss.pack(side='left')
lbl_dis_gauss = tk.Label(frm_gaussian, text="Distance")
lbl_dis_gauss.pack(side='left')
var_dis_gauss = tk.StringVar(root)  # variable for spinbox-value
var_dis_gauss.set(str(distance_gaussian))  # Initial value
spn_dis_gauss = tk.Spinbox(
    frm_gaussian, textvariable=var_dis_gauss, format="%.1f", from_=0., to=6.0, increment=0.1,
    command=lambda: change_distance_gaussian(float(var_dis_gauss.get())), width=5
)
spn_dis_gauss.pack(side='left')
btn_gaussian = tk.Button(frm_gaussian, text="Push", command=lambda: impact_gaussian())
btn_gaussian.pack(side='left')

# Play and pause button
btn_pp = tk.Button(root, text="Play/Pause", command=switch)
btn_pp.pack(side='left')

# Clear button
btn_clr = tk.Button(root, text="Reset", command=reset)
btn_clr.pack(side='left')

# Draw animation
set_axis()
anim = animation.FuncAnimation(fig, update, interval=50, save_count=100)
root.mainloop()
