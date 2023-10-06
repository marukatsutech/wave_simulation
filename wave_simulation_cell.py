# Wave simulation (Cell style)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
import matplotlib.ticker as ticker


def change_inc_dec_step(value):
    global step_inc_dec
    step_inc_dec = float(value)


def change_k(value):
    global k
    k = float(value)


def change_mass(value):
    global mass
    mass = float(value)


def clear_boundary():
    global boundary
    boundary = np.zeros((num_of_mass_x, num_of_mass_y))
    update_cells()


def clear_cells():
    global is_play, cnt, z, a, v, tx_step
    is_play = False
    cnt = 0
    z = np.zeros((num_of_mass_x, num_of_mass_y))
    a = np.zeros((num_of_mass_x, num_of_mass_y))
    v = np.zeros((num_of_mass_x, num_of_mass_y))
    update_cells()
    tx_step.set_text("Step=" + str(cnt))


def next_generation():
    global z, v
    for i in range(num_of_mass_y):
        for j in range(num_of_mass_x):
            v[j][i] = v[j][i] + a[j][i]
            z[j][i] = z[j][i] + v[j][i]


def z_boundary(j, i, j_neighbor, i_neighbor):
    if boundary[j_neighbor][i_neighbor] == 0:
        return z[j_neighbor][i_neighbor]
    else:
        if var_boundary_cond.get() == 1:  # Fixed end
            return 0.
        else:
            return z[j][i]


def eval_cells():
    global z, a
    for i in range(num_of_mass_y):
        for j in range(num_of_mass_x):
            if boundary[j][i] == 0:
                if var_boundary_cond.get() == 1:  # Fixed end
                    z_boundary_left = 0.
                    z_boundary_right = 0.
                    z_boundary_top = 0.
                    z_boundary_bottom = 0.
                else:  # Free end
                    z_boundary_left = z[0][i]
                    z_boundary_right = z[num_of_mass_x - 1][i]
                    z_boundary_top = z[j][num_of_mass_y - 1]
                    z_boundary_bottom = z[j][0]
                if j == 0:
                    diff_left = z[j][i] - z_boundary_left
                else:
                    diff_left = z[j][i] - z_boundary(j, i, j - 1, i)
                if j == num_of_mass_x - 1:
                    diff_right = z[j][i] - z_boundary_right
                else:
                    diff_right = z[j][i] - z_boundary(j, i, j + 1, i)
                if i == 0:
                    diff_bottom = z[j][i] - z_boundary_bottom
                else:
                    diff_bottom = z[j][i] - z_boundary(j, i, j, i - 1)
                if i == num_of_mass_y - 1:
                    diff_top = z[j][i] - z_boundary_top
                else:
                    diff_top = z[j][i] - z_boundary(j, i, j, i + 1)
                force = - k * (diff_left + diff_right + diff_bottom + diff_top)
                a[j][i] = force / mass


def update_cells():     # Update scat
    global scat, scat_boundary, x_cells0, y_cells0, s_cells0, c_cells0, x_boundary, y_boundary, s_boundary
    x_cells0.clear()
    y_cells0.clear()
    s_cells0.clear()
    c_cells0.clear()
    x_boundary.clear()
    y_boundary.clear()
    s_boundary.clear()
    for i in range(num_of_mass_y):
        for j in range(num_of_mass_x):
            if z[j][i] != 0.:
                y_cells0.append(i)
                x_cells0.append(j)
                s_cells0.append(size_maker)
                c_cells0.append(z[j][i])
            if boundary[j][i] == 1:
                y_boundary.append(i)
                x_boundary.append(j)
                s_boundary.append(size_maker)
    scat.set_offsets(np.column_stack([x_cells0, y_cells0]))
    scat.set_sizes(s_cells0)
    scat.set_array(np.array(c_cells0))
    scat_boundary.set_offsets(np.column_stack([x_boundary, y_boundary]))
    scat_boundary.set_sizes(s_boundary)


def mouse_motion(event):
    global lbl_value_cell
    if event.dblclick == 1:
        # print("double click")
        pass
    elif event.button == 1:
        # print("left click")
        if str(event.xdata) != "None" and str(event.ydata) != "None":
            # print(event.xdata, event.ydata)
            if 0 <= round(event.xdata) <= num_of_mass_x - 1 and 0 <= round(event.ydata) <= num_of_mass_y - 1:
                # print(round(event.xdata), round(event.ydata))
                if var_radio_click_option.get() == 1:
                    z[round(event.xdata)][round(event.ydata)] += step_inc_dec
                elif var_radio_click_option.get() == 2:
                    z[round(event.xdata)][round(event.ydata)] -= step_inc_dec
                elif var_radio_click_option.get() == 3:
                    boundary[round(event.xdata)][round(event.ydata)] = 1
                elif var_radio_click_option.get() == 4:
                    boundary[round(event.xdata)][round(event.ydata)] = 0
                else:
                    pass
                lbl_value_cell['text'] = " Cell(x" + str(round(event.xdata)) + ", y" + str(round(event.ydata)) \
                                         + ")=" + str(z[round(event.xdata)][round(event.ydata)])
                update_cells()
                eval_cells()
    elif event.button == 3:
        # print("right click")
        pass


def on_change_window(e):
    if not is_play:
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
x_min = - 1.
x_max = 50.
y_min = - 1.
y_max = 50.

# Animation control
cnt = 0
is_play = False

# Parameter
num_of_mass_x = 50
num_of_mass_y = 50
mass_init = 5.
mass = mass_init
k_init = 1.
k = k_init

step_inc_dec = 5

z = np.zeros((num_of_mass_x, num_of_mass_y))
a = np.zeros((num_of_mass_x, num_of_mass_y))
v = np.zeros((num_of_mass_x, num_of_mass_y))
boundary = np.zeros((num_of_mass_x, num_of_mass_y))

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111)
ax.set_title("Wave simulation (Cell style)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_aspect("equal")
ax.grid()
# ax.invert_yaxis()
x_ticks = []
for ix in range(num_of_mass_x):
    x_ticks.append(ix)
ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
y_ticks = []
for iy in range(num_of_mass_y):
    y_ticks.append(iy)
ax.yaxis.set_major_locator(ticker.FixedLocator(y_ticks))
ax.set_xticklabels(x_ticks, fontsize=6)
ax.set_yticklabels(y_ticks, fontsize=6)

# Generate graphic items
tx_step = ax.text(x_min, y_max * 0.95, "Step=" + str(0))
x_cells0 = []
y_cells0 = []
s_cells0 = []  # maker size
c_cells0 = []  # maker color
size_maker = 80
scat = ax.scatter(x_cells0, y_cells0, marker='s', s=size_maker, c=c_cells0, cmap='seismic', vmin=-5., vmax=5.)
x_boundary = []
y_boundary = []
s_boundary = []  # maker size
c_boundary = []  # maker color
size_maker = 80
scat_boundary = ax.scatter(x_boundary, y_boundary, marker='s', s=size_maker, c='black')


update_cells()

# Tkinter
root = tk.Tk()
root.title("Wave simulation (cell style)")
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
btn_clr = tk.Button(frm_ac, text="Clear cells", command=clear_cells)
btn_clr.pack()
# Clear boundary
btn_bd = tk.Button(frm_ac, text="Clear boundary", command=clear_boundary)
btn_bd.pack()

# Boundary condition
frm_bd = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
frm_bd.pack(side='left', fill=tk.Y)
var_boundary_cond = tk.IntVar(root)
rdb_fxe = tk.Radiobutton(frm_bd, text="Fixed end", value=1, var=var_boundary_cond)
rdb_fxe.pack()
rdb_fre = tk.Radiobutton(frm_bd, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack()
var_boundary_cond.set(1)

# Cells
frm_parameter = ttk.Labelframe(root, relief="ridge", text="Cells parameters", labelanchor="n")
frm_parameter.pack(side='left', fill=tk.Y)
lbl_mass = tk.Label(frm_parameter, text="Mass:")
lbl_mass.pack()
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_mass = tk.Spinbox(
    frm_parameter, textvariable=var_mass, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_mass(var_mass.get()), width=6
    )
spn_mass.pack()
lbl_k = tk.Label(frm_parameter, text="k(Constant of springs between cells):")
lbl_k.pack()
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    frm_parameter, textvariable=var_k, format="%.3f", from_=0.1, to=100., increment=0.01,
    command=lambda: change_k(var_k.get()), width=6
    )
spn_k.pack()

# Radio button
frm_cp = ttk.Labelframe(root, relief="ridge", text="Click option", labelanchor="n")
frm_cp.pack(side='left', fill=tk.Y)
var_radio_click_option = tk.IntVar(root)
# Radio button 1st
r_increase = tk.Radiobutton(frm_cp, text="Check value", value=0, var=var_radio_click_option)
r_increase.pack()
# Radio button 2st
r_increase = tk.Radiobutton(frm_cp, text="Increase", value=1, var=var_radio_click_option)
r_increase.pack()
# Radio button 3nd
r_decrease = tk.Radiobutton(frm_cp, text="Decrease", value=2, var=var_radio_click_option)
r_decrease.pack()
# Radio button 4th
r_boundary_add = tk.Radiobutton(frm_cp, text="Add boundary", value=3, var=var_radio_click_option)
r_boundary_add.pack()
# Radio button 5th
r_boundary_del = tk.Radiobutton(frm_cp, text="Delete boundary", value=4, var=var_radio_click_option)
r_boundary_del.pack()
var_radio_click_option.set(0)  # set default
# Label inc/dec step
lbl_step = tk.Label(frm_cp, text="Increase/decrease step:")
lbl_step.pack()
var_step = tk.StringVar(root)  # variable for spinbox-value
var_step.set(step_inc_dec)  # Initial value
spn_step = tk.Spinbox(
    frm_cp, textvariable=var_step, format="%.1f", from_=1, to=10, increment=1,
    command=lambda: change_inc_dec_step(var_step.get()), width=6
    )
spn_step.pack(side='left')

lbl_value_cell = tk.Label(root, text=" Cell(x0, y0)=0")
lbl_value_cell.pack(side='left')

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=100)
root.bind('<Configure>', on_change_window)
root.mainloop()
