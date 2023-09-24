# Wave simulation (Detail)(Displacement, velocity and acceleration)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk


def change_limit_diff(value):
    global limit_diff
    limit_diff = float(value)


def change_constraint_area_x(value):
    global constraint_location
    constraint_location = float(value)
    line_constraint_upper.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                     [constraint_height, constraint_height])
    line_constraint_lower.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                   [-constraint_height, -constraint_height])


def change_constraint_area_width(value):
    global constraint_width
    constraint_width = float(value)
    line_constraint_upper.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                     [constraint_height, constraint_height])
    line_constraint_lower.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                   [-constraint_height, -constraint_height])


def change_constraint_area_height(value):
    global constraint_height
    constraint_height = float(value)
    line_constraint_upper.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                     [constraint_height, constraint_height])
    line_constraint_lower.set_data([constraint_location - constraint_width, constraint_location + constraint_width],
                                   [-constraint_height, -constraint_height])


def change_range_x(value):
    global x_min, x_max, num_of_mass
    x_min = 0.
    x_max = float(value)
    num_of_mass = int(x_max)
    ax.set_xlim(x_min, x_max)


def change_range_y(value):
    global tx_step, y_min, y_max
    y_min = - float(value)
    y_max = float(value)
    ax.set_ylim(y_min, y_max)
    tx_step.set_position([0, y_max * 0.95])


def change_sigma(value):
    global sigma
    sigma = float(value)


def change_k1(value):
    global k1
    k1 = float(value)


def change_k1_op(value):
    global k1_op
    k1_op = float(value)


def change_k0_op(value):
    global k0_op
    k0_op = float(value)


def change_k(value):
    global k
    k = float(value)


def change_mass(value):
    global mass
    mass = float(value)


def clear_cells():
    global is_play, cnt, tx_step
    is_play = False
    cnt = 0
    tx_step.set_text("Step=" + str(cnt))
    for i in range(num_of_mass_max):
        y[i] = 0.
        dl[i] = 0.
        dr[i] = 0.
        f[i] = 0.
        a[i] = 0.
        v[i] = 0.
    update_cells()


def next_generation():
    global y, v
    for i in range(num_of_mass):
        v[i] = v[i] + a[i]
        y[i] = y[i] + v[i]


def eval_cells():
    global y, dl, dr, f, a, v
    if var_boundary_cond.get() == 1:    # Fixed end
        y_boundary_left = 0.
        y_boundary_right = 0.
    else:                               # Free end
        y_boundary_left = y[0]
        y_boundary_right = y[num_of_mass - 1]
    for i in range(num_of_mass):
        k0_current = k
        if var_rd_k0_op.get():
            if constraint_location - constraint_width <= i <= constraint_location + constraint_width:
                k0_current = k0_op
            else:
                k0_current = k
        if i == 0:
            dl[i] = y[i] - y_boundary_left
        else:
            dl[i] = y[i] - y[i - 1]
        if i == num_of_mass - 1:
            dr[i] = y[i] - y_boundary_right
        else:
            dr[i] = y[i] - y[i + 1]
        if var_rd_k1_op.get() == 0:
            f[i] = - k0_current * (dl[i] + dr[i])
        elif var_rd_k1_op.get() == 1:
            f[i] = - k0_current * (dl[i] + dr[i]) - k1 * y[i]
        elif var_rd_k1_op.get() == 2:
            if constraint_location - constraint_width <= i <= constraint_location + constraint_width:
                f[i] = - k0_current * (dl[i] + dr[i]) - k1 * y[i]
            else:
                f[i] = - k0_current * (dl[i] + dr[i])
        elif var_rd_k1_op.get() == 3:
            if constraint_location - constraint_width <= i <= constraint_location + constraint_width:
                if - constraint_height <= y[i] <= constraint_height:
                    f[i] = - k0_current * (dl[i] + dr[i])
                elif y[i] > constraint_height:
                    f[i] = - k0_current * (dl[i] + dr[i]) - k1 * (y[i] - constraint_height)
                else:
                    f[i] = - k0_current * (dl[i] + dr[i]) - k1 * (y[i] + constraint_height)
            else:
                f[i] = - k0_current * (dl[i] + dr[i])
        elif var_rd_k1_op.get() == 4:
            if constraint_location - constraint_width <= i <= constraint_location + constraint_width:
                over_f_l = 0.
                over_f_r = 0.
                if np.abs(dl[i]) > limit_diff:
                    over_f_l = - k1 * np.sign(dl[i]) * (np.abs(dl[i]) - limit_diff)
                if np.abs(dr[i]) > limit_diff:
                    over_f_r = - k1 * np.sign(dr[i]) * (np.abs(dr[i]) - limit_diff)
                f[i] = - k0_current * (dl[i] + dr[i]) + over_f_l + over_f_r
            else:
                f[i] = - k0_current * (dl[i] + dr[i])
        else:
            if constraint_location - constraint_width <= i <= constraint_location + constraint_width:
                k1_current = k1_op
            else:
                k1_current = k1
            f[i] = - k0_current * (dl[i] + dr[i]) - k1_current * y[i]
        a[i] = f[i] / mass


def update_cells():
    global string_y, string_v, string_a
    # draw string
    string_y.set_ydata(y)
    string_v.set_ydata(v)
    string_a.set_ydata(a)


def impact_gauss(impact_x, impact_y):
    global y, v, a
    if var_round.get():
        h_gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- 0. ** 2 / (2 * sigma ** 2))
        gaussian = np.round(np.round(1 / np.sqrt(2 * np.pi * sigma ** 2) *
                                     np.e ** (- (x - round(impact_x)) ** 2 / (2 * sigma ** 2)) * impact_y / h_gaussian))
    else:
        h_gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- 0. ** 2 / (2 * sigma ** 2))
        gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- (x - impact_x) ** 2 / (2 * sigma ** 2)) * \
               impact_y / h_gaussian
    if var_radio_y_or_v.get() == 0:
        y = y.copy() + gaussian
    elif var_radio_y_or_v.get() == 1:
        v = v.copy() + gaussian
    else:
        pass
    if var_boundary_cond.get() == 1:
        y[0] = 0.
        y[len(x) - 1] = 0.
        v[0] = 0.
        v[len(x) - 1] = 0.
        a[0] = 0.
        a[len(x) - 1] = 0.


def mouse_motion(event):
    if event.dblclick == 1:
        # print("double click")
        pass
    elif event.button == 1:
        # print("left click")
        if str(event.xdata) != "None" and str(event.ydata) != "None":
            # print(event.xdata, event.ydata)
            if 0 <= round(event.xdata) <= num_of_mass - 1:
                impact_gauss(event.xdata, event.ydata)
                update_cells()
                eval_cells()
    elif event.button == 3:
        # print("right click")
        pass


def step():
    global cnt
    cnt += 1
    tx_step.set_text("Step=" + str(cnt))
    next_generation()
    eval_cells()
    update_cells()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True
    tx_step.set_text("Step=" + str(cnt))


def update(f):
    global tx_step, cnt
    if is_play:
        tx_step.set_text("Step=" + str(cnt))
        cnt += 1
        next_generation()
        eval_cells()
        update_cells()


# Global variables
range_x_init = 400.
range_y_init = 10.
x_min = 0.
x_max = range_x_init

y_min = - range_y_init
y_max = range_y_init

# Animation control
cnt = 0
is_play = False

# Parameter
num_of_mass_max = 800
num_of_mass = 400

mass_init = 10.
mass = mass_init
k_init = 1.
k = k_init
k0_op_init = 1.
k0_op = k0_op_init
k1_init = 1.
k1 = k1_init
k0_op_init = 1.
k0_op = k0_op_init
k1_op_init = 1.
k1_op = k0_op_init

sigma_init = 4.0
sigma = sigma_init

x = np.arange(num_of_mass_max)
y = x * 0.
dl = x * 0.
dr = x * 0.
f = x * 0.
a = x * 0.
v = x * 0.

constraint_location = (x_max - x_min) / 2.
constraint_width = 10.
constraint_height = 1.
limit_diff_init = 0.01
limit_diff = limit_diff_init

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111)
ax.set_title("Wave simulation (Displacement, velocity and acceleration)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)


# ax.set_aspect("equal")
ax.grid()
# ax.invert_yaxis()

# Generate graphic items
tx_step = ax.text(x_min, y_max * 0.95, "Step=" + str(0))
string_y, = ax.plot(x, y, linestyle='-', label="Displacement")
string_v, = ax.plot(x, y, linestyle='--', label="velocity", linewidth=1)
string_a, = ax.plot(x, y, linestyle='-.', label="acceleration", linewidth=1)
ax.legend(prop={"size": 8}, loc="best")

line_constraint_upper, = ax.plot([constraint_location - constraint_width, constraint_location + constraint_width],
                                 [constraint_height, constraint_height], color='red', linestyle='--')
line_constraint_lower, = ax.plot([constraint_location - constraint_width, constraint_location + constraint_width],
                                 [-constraint_height, -constraint_height], color='red', linestyle='--')

# Tkinter
root = tk.Tk()
root.title("Wave simulation (Digitized detail)")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')
canvas.mpl_connect('button_press_event', mouse_motion)

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

# x, y
frm_range = ttk.Labelframe(root, relief="ridge", text="Range", labelanchor="n", width=100)
frm_range.pack(side='left', fill=tk.Y)
lbl_x = tk.Label(frm_range, text="x(number of mass):")
lbl_x.pack()
var_x = tk.StringVar(root)  # variable for spinbox-value
var_x.set(range_x_init)  # Initial value
spn_x = tk.Spinbox(
    frm_range, textvariable=var_x, format="%.1f", from_=50, to=800, increment=50,
    command=lambda: change_range_x(var_x.get()), width=4
    )
spn_x.pack()
lbl_y = tk.Label(frm_range, text="y:")
lbl_y.pack()
var_y = tk.StringVar(root)  # variable for spinbox-value
var_y.set(range_y_init)  # Initial value
spn_y = tk.Spinbox(
    frm_range, textvariable=var_y, format="%.1f", from_=1., to=20.0, increment=1.,
    command=lambda: change_range_y(var_y.get()), width=4
    )
spn_y.pack()

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
lbl_mass = tk.Label(frm_parameter, text="Mass:")
lbl_mass.pack()
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_mass = tk.Spinbox(
    frm_parameter, textvariable=var_mass, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_mass(var_mass.get()), width=5
    )
spn_mass.pack()
lbl_k = tk.Label(frm_parameter, text="k0(between cells):")
lbl_k.pack()
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    frm_parameter, textvariable=var_k, format="%.1f", from_=0., to=100., increment=1.,
    command=lambda: change_k(var_k.get()), width=5
    )
spn_k.pack()
lbl_k1 = tk.Label(frm_parameter, text="k1(option):")
lbl_k1.pack()
var_k1 = tk.StringVar(root)  # variable for spinbox-value
var_k1.set(k1_init)  # Initial value
spn_k1 = tk.Spinbox(
    frm_parameter, textvariable=var_k1, format="%.1f", from_=0., to=100., increment=1.,
    command=lambda: change_k1(var_k1.get()), width=5
    )
spn_k1.pack()

# k0 option
frm_k0_op = ttk.Labelframe(root, relief="ridge", text="k0 option", labelanchor="n")
frm_k0_op.pack(side='left', fill=tk.Y)
var_rd_k0_op = tk.IntVar(root)
# Radio button 1st
rd_k0_op_non = tk.Radiobutton(frm_k0_op, text="none", value=0, var=var_rd_k0_op)
rd_k0_op_non.pack()
# Radio button 2st
rd_k0_op_all = tk.Radiobutton(frm_k0_op, text="Apply in constraint area", value=1, var=var_rd_k0_op)
rd_k0_op_all.pack()
var_rd_k0_op.set(0)
lbl_k0_op_ca2 = tk.Label(frm_k0_op, text="k0 in constraint area:")
lbl_k0_op_ca2.pack()
var_k0_op = tk.StringVar(root)  # variable for spinbox-value
var_k0_op.set(k0_op_init)  # Initial value
spn_k0_op = tk.Spinbox(
    frm_k0_op, textvariable=var_k0_op, format="%.2f", from_=0., to=10., increment=0.01,
    command=lambda: change_k0_op(var_k0_op.get()), width=5
    )
spn_k0_op.pack()

# k1 option
frm_k1_op = ttk.Labelframe(root, relief="ridge", text="k1 option", labelanchor="n")
frm_k1_op.pack(side='left', fill=tk.Y)
var_rd_k1_op = tk.IntVar(root)
# Radio button 1st
rd_k1_op_non = tk.Radiobutton(frm_k1_op, text="none", value=0, var=var_rd_k1_op)
rd_k1_op_non.pack()
# Radio button 2st
rd_k1_op_all = tk.Radiobutton(frm_k1_op, text="All area", value=1, var=var_rd_k1_op)
rd_k1_op_all.pack()
# Radio button 3nd
rd_k1_op_ca = tk.Radiobutton(frm_k1_op, text="Constraint area", value=2, var=var_rd_k1_op)
rd_k1_op_ca.pack()
# Radio button 4th
rd_k1_op_ca1 = tk.Radiobutton(frm_k1_op, text="Constraint area(over height limit)", value=3, var=var_rd_k1_op)
rd_k1_op_ca1.pack()
# Radio button 5th
rd_k1_op_ca2 = tk.Radiobutton(frm_k1_op, text="Constraint area(over limited diff.)", value=4, var=var_rd_k1_op)
rd_k1_op_ca2.pack()
lbl_k1_op_ca2 = tk.Label(frm_k1_op, text="Limit diff. between cells:")
lbl_k1_op_ca2.pack()
var_k1_diff = tk.StringVar(root)  # variable for spinbox-value
var_k1_diff.set(limit_diff_init)  # Initial value
spn_k1_diff = tk.Spinbox(
    frm_k1_op, textvariable=var_k1_diff, format="%.2f", from_=0., to=10., increment=0.01,
    command=lambda: change_limit_diff(var_k1_diff.get()), width=5
    )
spn_k1_diff.pack()
# Radio button 6th
rd_k1_op_ca2 = tk.Radiobutton(frm_k1_op, text="All area and apply following", value=5, var=var_rd_k1_op)
rd_k1_op_ca2.pack()
var_rd_k1_op.set(0)
lbl_k1_op_ca3 = tk.Label(frm_k1_op, text="k1 in constraint area:")
lbl_k1_op_ca3.pack()
var_k1_op = tk.StringVar(root)  # variable for spinbox-value
var_k1_op.set(k1_op_init)  # Initial value
spn_k1_op = tk.Spinbox(
    frm_k1_op, textvariable=var_k1_op, format="%.2f", from_=0., to=10., increment=0.01,
    command=lambda: change_k1_op(var_k1_op.get()), width=5
    )
spn_k1_op.pack()

# Constraint area
frm_const = ttk.Labelframe(root, relief="ridge", text="Constraint area", labelanchor="n")
frm_const.pack(side='left', fill=tk.Y)
lbl_x_area = tk.Label(frm_const, text="x:")
lbl_x_area.pack()
var_x_area = tk.StringVar(root)  # variable for spinbox-value
var_x_area.set(constraint_location)  # Initial value
spn_x_area = tk.Spinbox(
    frm_const, textvariable=var_x_area, format="%.2f", from_=0., to=800., increment=1.,
    command=lambda: change_constraint_area_x(var_x_area.get()), width=5
    )
spn_x_area.pack()
lbl_width_area = tk.Label(frm_const, text="Width:")
lbl_width_area.pack()
var_width_area = tk.StringVar(root)  # variable for spinbox-value
var_width_area.set(constraint_width)  # Initial value
spn_width_area = tk.Spinbox(
    frm_const, textvariable=var_width_area, format="%.2f", from_=1., to=400., increment=1.,
    command=lambda: change_constraint_area_width(var_width_area.get()), width=5
    )
spn_width_area.pack()
lbl_height_area = tk.Label(frm_const, text="Height")
lbl_height_area.pack()
var_height_area = tk.StringVar(root)  # variable for spinbox-value
var_height_area.set(constraint_height)  # Initial value
spn_height_area = tk.Spinbox(
    frm_const, textvariable=var_height_area, format="%.2f", from_=0., to=10., increment=0.1,
    command=lambda: change_constraint_area_height(var_height_area.get()), width=5
    )
spn_height_area.pack()

# Click option of Gaussian displacement
frm_click = ttk.Labelframe(root, relief="ridge", text="Click option(apply Gaussian)", labelanchor="n")
frm_click.pack(side='left', fill=tk.Y)
var_radio_y_or_v = tk.IntVar(root)
# Radio button 1st
r_y = tk.Radiobutton(frm_click, text="y(Displacement)", value=0, var=var_radio_y_or_v)
r_y.pack()
# Radio button 2nd
r_v = tk.Radiobutton(frm_click, text="v(Velocity)", value=1, var=var_radio_y_or_v)
r_v.pack()

# Gaussian
frm_gauss = ttk.Labelframe(root, relief="ridge", text="Gaussian", labelanchor="n")
frm_gauss.pack(side='left', fill=tk.Y)
lbl_sigma = tk.Label(frm_gauss, text="Sigma:")
lbl_sigma.pack()
var_sgm = tk.StringVar(root)  # variable for spinbox-value
var_sgm.set(sigma_init)  # Initial value
spn_sgm = tk.Spinbox(
    frm_gauss, textvariable=var_sgm, format="%.1f", from_=0.1, to=10.0, increment=0.1,
    command=lambda: change_sigma(var_sgm.get()), width=4
    )
spn_sgm.pack()
# Checkbutton of round on/off
var_round = tk.BooleanVar(root)    # Variable for checkbutton
chk_round = tk.Checkbutton(frm_gauss, text="Round", variable=var_round)
chk_round.pack()

update_cells()

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=100)
root.mainloop()
