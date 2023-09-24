# Wave simulation (Spiral)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d
import tkinter as tk
from tkinter import ttk


def change_range_x(value):
    global x_min, x_max, num_of_moment
    x_min = 0.
    x_max = float(value)
    num_of_moment = int(x_max)
    ax.set_xlim(x_min, x_max)


def change_range_y(value):
    global tx_step, y_min, y_max
    y_min = - float(value)
    y_max = float(value)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(y_min, y_max)
    tx_step.set_position([0, y_max, z_max * 1.2])


def change_loc_gauss(self):
    global loc_gaussian, line_loc, data_line_loc
    loc_gaussian = x_min + (x_max - x_min) * var_scale.get() / 100.
    # print(self, loc_gaussian)
    data_line_loc = np.array([[loc_gaussian, loc_gaussian], [0., 0.], [0., z_min]])
    line_loc.set_data(data_line_loc[0, 0:2], data_line_loc[1, 0:2])
    line_loc.set_3d_properties(data_line_loc[2, 0:2])


def change_modulo(value):
    global modulo
    modulo = float(value)


def change_rot_angle_v(value):
    global rot_angle_v
    rot_angle_v = float(value)
    
    
def change_angle_gauss(value):
    global angle_gauss
    angle_gauss = float(value)


def change_sigma(value):
    global sigma
    sigma = float(value)


def change_k1(value):
    global k1
    k1 = float(value)


def change_k(value):
    global k
    k = float(value)


def change_moment(value):
    global moment_i
    moment_i = float(value)


def clear_cells():
    global is_play, cnt, tx_step, angle, dl, dr, torque, angle_a, angle_v, angle_input, is_rot
    is_play = False
    cnt = 0
    tx_step.set_text("Step=" + str(cnt))
    angle_input = 0.
    is_rot = False
    for i in range(num_of_moment_max):
        angle[i] = 0.
        dl[i] = 0.
        dr[i] = 0.
        torque[i] = 0.
        angle_a[i] = 0.
        angle_v[i] = 0.
    update_cells()


def next_generation():
    global angle, angle_v, angle_input
    for i in range(num_of_moment):
        angle_v[i] = angle_v[i] + angle_a[i]
        angle[i] = angle[i] + angle_v[i]
    if is_rot:
        angle_input += rot_angle_v * np.pi
        angle[int(loc_gaussian)] = angle_input


def eval_cells_modulo():
    global angle, dl, dr, torque, angle_a, angle_v
    if var_boundary_cond.get() == 1:    # Fixed end
        angle_boundary_left = 0.
        angle_boundary_right = 0.
    else:                               # Free end
        angle_boundary_left = angle[0]
        angle_boundary_right = angle[num_of_moment - 1]
    for i in range(num_of_moment):
        if i == 0:
            dl[i] = np.sign(angle[i]) * (np.abs(angle[i]) % modulo) - np.sign(angle_boundary_left) * \
                    (np.abs(angle_boundary_left) % modulo)
        else:
            dl[i] = np.sign(angle[i]) * (np.abs(angle[i]) % modulo) - np.sign(angle[i - 1]) * \
                    (np.abs(angle[i - 1]) % modulo)
        if i == num_of_moment - 1:
            dr[i] = np.sign(angle[i]) * (np.abs(angle[i]) % modulo) - np.sign(angle_boundary_right) * \
                    (np.abs(angle_boundary_right) % modulo)
        else:
            dr[i] = np.sign(angle[i]) * (np.abs(angle[i]) % modulo) - np.sign(angle[i + 1]) * \
                    (np.abs(angle[i + 1]) % modulo)
        if var_chk_k1.get():
            torque[i] = - k * (dl[i] + dr[i]) - k1 * angle[i]
        else:
            torque[i] = - k * (dl[i] + dr[i])
        angle_a[i] = torque[i] / moment_i


def eval_cells():
    global angle, dl, dr, torque, angle_a, angle_v
    if var_boundary_cond.get() == 1:    # Fixed end
        angle_boundary_left = 0.
        angle_boundary_right = 0.
    else:                               # Free end
        angle_boundary_left = angle[0]
        angle_boundary_right = angle[num_of_moment - 1]
    for i in range(num_of_moment):
        if i == 0:
            dl[i] = angle[i] - angle_boundary_left
        else:
            dl[i] = angle[i] - angle[i - 1]
        if i == num_of_moment - 1:
            dr[i] = angle[i] - angle_boundary_right
        else:
            dr[i] = angle[i] - angle[i + 1]
        if var_chk_k1.get():
            torque[i] = - k * (dl[i] + dr[i]) - k1 * angle[i]
        else:
            torque[i] = - k * (dl[i] + dr[i])
        angle_a[i] = torque[i] / moment_i


def update_cells():
    global string_angle, angle, y, z, string_angle_not_rolled, string_v, string_a
    # draw string
    z = np.sin(angle * np.pi)
    y = np.cos(angle * np.pi)
    data_angle = np.array([x, y, z])
    string_angle.set_data(data_angle[0, 0:num_of_moment], data_angle[1, 0:num_of_moment])
    string_angle.set_3d_properties(data_angle[2, 0:num_of_moment])
    data_z = np.array([x, x * 0. + y_max, z])
    string_z.set_data(data_z[0, 0:num_of_moment], data_z[1, 0:num_of_moment])
    string_z.set_3d_properties(data_z[2, 0:num_of_moment])
    data_angle_nr = np.array([x, angle, x * 0. + z_min])
    string_angle_not_rolled.set_data(data_angle_nr[0, 0:num_of_moment], data_angle_nr[1, 0:num_of_moment])
    string_angle_not_rolled.set_3d_properties(data_angle_nr[2, 0:num_of_moment])
    data_v = np.array([x, angle_v, x * 0. + z_min])
    string_v.set_data(data_v[0, 0:num_of_moment], data_v[1, 0:num_of_moment])
    string_v.set_3d_properties(data_v[2, 0:num_of_moment])
    data_a = np.array([x, angle_a, x * 0. + z_min])
    string_a.set_data(data_a[0, 0:num_of_moment], data_a[1, 0:num_of_moment])
    string_a.set_3d_properties(data_a[2, 0:num_of_moment])


def input_gauss():
    global angle, angle_v, angle_a
    h_gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- 0. ** 2 / (2 * sigma ** 2))
    gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- (x - loc_gaussian) ** 2 / (2 * sigma ** 2)) * \
               angle_gauss / h_gaussian
    angle = angle.copy() + gaussian
    if var_boundary_cond.get() == 1:
        angle[0] = 0.
        angle[len(x) - 1] = 0.
        angle_v[0] = 0.
        angle_v[len(x) - 1] = 0.
        angle_a[0] = 0.
        angle_a[len(x) - 1] = 0.
    update_cells()


def input_gauss_rot():
    global is_rot, is_play
    if var_input_op.get() == 0:
        input_gauss()
    else:
        if is_rot:
            is_rot = False
        else:
            is_rot = True
            is_play = True


def step():
    global cnt
    cnt += 1
    tx_step.set_text("Step=" + str(cnt))

    if var_modulo.get():
        eval_cells()
    else:
        eval_cells()
    next_generation()
    update_cells()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True
    tx_step.set_text("Step=" + str(cnt))


def update(f):
    global tx_step, cnt, angle
    if is_play:
        tx_step.set_text("Step=" + str(cnt))
        cnt += 1
        if var_modulo_op.get() == 0:
            eval_cells()
            next_generation()
        elif var_modulo_op.get() == 1:
            eval_cells()
            next_generation()
            angle = np.sign(angle) * (np.abs(angle) % modulo)
        else:
            eval_cells_modulo()
            next_generation()
        update_cells()


# Global variables
range_x_init = 400.
range_y_init = 1.5
range_z_init = 1.5

x_min = 0.
x_max = range_x_init
y_min = - range_y_init
y_max = range_y_init
z_min = - range_z_init
z_max = range_z_init

# Animation control
cnt = 0
is_play = False

# Parameter
num_of_moment_max = 800
num_of_moment = 400
moment_i_init = 2.
moment_i = moment_i_init
k_init = 1.
k = k_init
k1_init = 1.
k1 = k1_init

sigma_init = 6.0
sigma = sigma_init

rot_angle_v_init = 0.01
rot_angle_v = rot_angle_v_init
is_rot = False
angle_input = 0.

modulo_init = 2.
modulo = modulo_init

x = np.arange(num_of_moment_max)
y = x * 0. + 1.
z = x * 0.

angle = x * 0.

dl = x * 0.
dr = x * 0.
torque = x * 0.
angle_a = x * 0.
angle_v = x * 0.

loc_gaussian = (x_max - x_min) / 2.
angle_gauss = 1.

# Generate figure and axes
fig = Figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
# ax = p3.Axes3D(fig)
ax.set_box_aspect((2, 1, 1))

ax.set_title("Wave simulation (Spiral)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel('z')
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_zlim(z_min, z_max)


# ax.set_aspect("equal")
ax.grid()
# ax.invert_yaxis()

# Generate graphic items
tx_step = ax.text(0, y_max, z_max * 1.2, "Step=" + str(0))

string_angle, = ax.plot(x, y, z, linestyle='-', label="Angle")
string_z, = ax.plot(x, y, z, linestyle='-', color="gray", label="z component", linewidth=1)
string_angle_not_rolled, = ax.plot(x, y, z, linestyle='--', label="Angle(not rolled)", linewidth=1)
string_v, = ax.plot(x, y, z, linestyle='--', label="Angle velocity", linewidth=1)
string_a, = ax.plot(x, y, z, linestyle='-.', label="Angle acceleration", linewidth=1)
ax.legend(prop={"size": 8}, loc="best")

data_line_loc = np.array([[loc_gaussian, loc_gaussian], [0., 0.], [0., z_min]])
line_loc, = ax.plot(data_line_loc[0, 0:2], data_line_loc[1, 0:2], data_line_loc[2, 0:2],
                    linestyle='--', linewidth=1, color='red')

string_center, = ax.plot(x, x * 0., x * 0., linestyle='-.', color="gray", linewidth=1)
c = Circle((0, 0), 1, ec='gray', fill=False, linestyle='--')
ax.add_patch(c)
art3d.pathpatch_2d_to_3d(c, z=0, zdir="x")
c1 = Circle((0, 0), 1, ec='gray', fill=False, linestyle='--')
ax.add_patch(c1)
art3d.pathpatch_2d_to_3d(c1, z=x_max, zdir="x")

# Tkinter`
root = tk.Tk()
root.title("Wave simulation (Spiral)")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Animation control
frm_ac = ttk.Labelframe(root, relief="ridge", text="Animation control", labelanchor="n", width=100)
frm_ac.pack(side='left', fill=tk.Y)
# Play and pause button
btn_pp = tk.Button(frm_ac, text="Play/Pause", command=switch)
btn_pp.pack()
# Step button
btn_pp = tk.Button(frm_ac, text="Step", command=step)
btn_pp.pack()
# Clear button
btn_clr = tk.Button(frm_ac, text="Clear", command=clear_cells)
btn_clr.pack()

# Range of x
frm_range = ttk.Labelframe(root, relief="ridge", text="Range", labelanchor="n", width=100)
frm_range.pack(side='left', fill=tk.Y)
lbl_x = tk.Label(frm_range, text="x(number of mass):")
lbl_x.pack()
var_x = tk.StringVar(root)  # variable for spinbox-value
var_x.set(range_x_init)  # Initial value
spn_x = tk.Spinbox(
    frm_range, textvariable=var_x, format="%.1f", from_=50, to=800, increment=50,
    command=lambda: change_range_x(var_x.get()), width=6
    )
spn_x.pack()


# Boundary condition
frm_bd = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
frm_bd.pack(side='left', fill=tk.Y)
var_boundary_cond = tk.IntVar(root)
rdb_fxe = tk.Radiobutton(frm_bd, text="Fixed end", value=1, var=var_boundary_cond)
rdb_fxe.pack()
rdb_fre = tk.Radiobutton(frm_bd, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack()
var_boundary_cond.set(1)

# Cells parameters
frm_parameter = ttk.Labelframe(root, relief="ridge", text="Cells parameters", labelanchor="n")
frm_parameter.pack(side='left', fill=tk.Y)
lbl_moment = tk.Label(frm_parameter, text="moment of inertia:")
lbl_moment.pack()
var_moment = tk.StringVar(root)  # variable for spinbox-value
var_moment.set(moment_i_init)  # Initial value
spn_moment = tk.Spinbox(
    frm_parameter, textvariable=var_moment, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_moment(var_moment.get()), width=6
    )
spn_moment.pack()
lbl_k = tk.Label(frm_parameter, text="k0(between cells):")
lbl_k.pack()
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    frm_parameter, textvariable=var_k, format="%.1f", from_=0., to=100., increment=1.,
    command=lambda: change_k(var_k.get()), width=6
    )
spn_k.pack()
# Checkbutton of K1 on/off
var_chk_k1 = tk.BooleanVar(root)    # Variable for checkbutton
chk_k1 = tk.Checkbutton(frm_parameter, text="k1(against angle):", variable=var_chk_k1)
chk_k1.pack()
var_k1 = tk.StringVar(root)  # variable for spinbox-value
var_k1.set(k1_init)  # Initial value
spn_k1 = tk.Spinbox(
    frm_parameter, textvariable=var_k1, format="%.2f", from_=0., to=100., increment=0.01,
    command=lambda: change_k1(var_k1.get()), width=6
    )
spn_k1.pack()

# Modulo
frm_modulo = ttk.Labelframe(root, relief="ridge", text="Modulo between cells", labelanchor="n")
frm_modulo.pack(side='left', fill=tk.Y)
var_modulo_op = tk.IntVar(root)
rdb_mod_none = tk.Radiobutton(frm_modulo, text="None", value=0, var=var_modulo_op)
rdb_mod_none.pack()
rdb_mod_angle = tk.Radiobutton(frm_modulo, text="Angle", value=1, var=var_modulo_op)
rdb_mod_angle.pack()
rdb_mod_diff = tk.Radiobutton(frm_modulo, text="Angle between cells", value=2, var=var_modulo_op)
rdb_mod_diff.pack()
lbl_modulo = tk.Label(frm_modulo, text="Quotient:")
lbl_modulo.pack()
var_modulo_op.set(0)
var_modulo = tk.StringVar(root)  # variable for spinbox-value
var_modulo.set(modulo_init)  # Initial value
spn_modulo = tk.Spinbox(
    frm_modulo, textvariable=var_modulo, format="%.1f", from_=0.1, to=100.0, increment=0.1,
    command=lambda: change_modulo(var_modulo.get()), width=6
    )
spn_modulo.pack()

# Input option
frm_input = ttk.Labelframe(root, relief="ridge", text="Input option", labelanchor="n")
frm_input.pack(side='left', fill=tk.Y)
var_input_op = tk.IntVar(root)
rdb_input_gauss = tk.Radiobutton(frm_input, text="Gaussian", value=0, var=var_input_op)
rdb_input_gauss.pack()
rdb_input_rot = tk.Radiobutton(frm_input, text="Rotation", value=1, var=var_input_op)
rdb_input_rot.pack()
var_input_op.set(0)
lbl_position = tk.Label(frm_input, text="Position:")
lbl_position.pack()
var_scale = tk.DoubleVar(root)
scl = tk.Scale(frm_input, variable=var_scale, orient='horizontal', length=100, command=change_loc_gauss)
scl.pack()
var_scale.set(50)
btn_input = tk.Button(frm_input, text="Input/start", command=input_gauss_rot)
btn_input .pack()

# Gaussian curve
frm_gauss = ttk.Labelframe(root, relief="ridge", text="Gaussian curve", labelanchor="n")
frm_gauss.pack(side='left', fill=tk.Y)
lbl_sigma = tk.Label(frm_gauss, text="sigma:")
lbl_sigma.pack()
var_sgm = tk.StringVar(root)  # variable for spinbox-value
var_sgm.set(sigma_init)  # Initial value
spn_sgm = tk.Spinbox(
    frm_gauss, textvariable=var_sgm, format="%.1f", from_=0.1, to=50.0, increment=0.1,
    command=lambda: change_sigma(var_sgm.get()), width=6
    )
spn_sgm.pack()
lbl_ang = tk.Label(frm_gauss, text="Angle(pi):")
lbl_ang.pack()
var_ang = tk.StringVar(root)  # variable for spinbox-value
var_ang.set(angle_gauss)  # Initial value
spn_ang = tk.Spinbox(
    frm_gauss, textvariable=var_ang, format="%.1f", from_=0.1, to=10.0, increment=0.1,
    command=lambda: change_angle_gauss(var_ang.get()), width=6
    )
spn_ang.pack()

# Rotation
frm_rot = ttk.Labelframe(root, relief="ridge", text="Rotation", labelanchor="n")
frm_rot.pack(side='left', fill=tk.Y)
lbl_ang_v = tk.Label(frm_rot, text="Angle velocity(pi):")
lbl_ang_v.pack()
var_rot = tk.StringVar(root)  # variable for spinbox-value
var_rot.set(rot_angle_v_init)  # Initial value
spn_rot = tk.Spinbox(
    frm_rot, textvariable=var_rot, format="%.3f", from_=0., to=50.0, increment=0.001,
    command=lambda: change_rot_angle_v(var_rot.get()), width=6
    )
spn_rot.pack()

update_cells()

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=100)
root.mainloop()
