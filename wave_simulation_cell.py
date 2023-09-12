# Wave simulation (Cell)

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


def eval_cells():
    global z, a
    for i in range(num_of_mass_y):
        for j in range(num_of_mass_x):
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
                diff_left = z[j][i] - z[j - 1][i]
            if j == num_of_mass_x - 1:
                diff_right = z[j][i] - z_boundary_right
            else:
                diff_right = z[j][i] - z[j + 1][i]
            if i == 0:
                diff_bottom = z[j][i] - z_boundary_bottom
            else:
                diff_bottom = z[j][i] - z[j][i - 1]
            if i == num_of_mass_y - 1:
                diff_top = z[j][i] - z_boundary_top
            else:
                diff_top = z[j][i] - z[j][i + 1]
            force = - k * (diff_left + diff_right + diff_bottom + diff_top)
            a[j][i] = force / mass


def update_cells():     # Update scat
    global x_cells0, y_cells0, s_cells0, c_cells0
    x_cells0.clear()
    y_cells0.clear()
    s_cells0.clear()
    c_cells0.clear()
    for i in range(num_of_mass_y):
        for j in range(num_of_mass_x):
            if z[j][i] != 0.:
                y_cells0.append(i)
                x_cells0.append(j)
                s_cells0.append(size_maker)
                c_cells0.append(z[j][i])
    scat.set_offsets(np.column_stack([x_cells0, y_cells0]))
    scat.set_sizes(s_cells0)
    scat.set_array(np.array(c_cells0))


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
                if var_radio_increase_decrease.get() == 1:
                    z[round(event.xdata)][round(event.ydata)] += step_inc_dec
                elif var_radio_increase_decrease.get() == 2:
                    z[round(event.xdata)][round(event.ydata)] -= step_inc_dec
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
mass = 1.
mass_init = mass
k = 0.5
k_init = k

step_inc_dec = 1

z = np.zeros((num_of_mass_x, num_of_mass_y))
a = np.zeros((num_of_mass_x, num_of_mass_y))
v = np.zeros((num_of_mass_x, num_of_mass_y))

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111)
ax.set_title("Wave simulation (Cell)")
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
x_cells0 = []
y_cells0 = []
s_cells0 = []  # maker size
c_cells0 = []  # maker color
size_maker = 80
scat = ax.scatter(x_cells0, y_cells0, marker='s', s=size_maker, c=c_cells0, cmap='seismic', vmin=-5., vmax=5.)
tx_step = ax.text(x_min, y_max * 0.95, "Step=" + str(0))

update_cells()

# Tkinter
root = tk.Tk()
root.title("Wave simulation (cell)")
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill='both')
canvas.mpl_connect('button_press_event', mouse_motion)

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

# Play and pause button
btn_pp = tk.Button(root, text="Play/Pause", command=switch)
btn_pp.pack(side='left')

# Step button
btn_pp = tk.Button(root, text="Step", command=step)
btn_pp.pack(side='left')

# Clear button
btn_clr = tk.Button(root, text="Clear", command=clear_cells)
btn_clr.pack(side='left')

# Boundary condition
frm1 = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
frm1.pack(side='left')
var_boundary_cond = tk.IntVar(root)
rdb_fxe = tk.Radiobutton(frm1, text="Fixed end, ", value=1, var=var_boundary_cond)
rdb_fxe.pack(side='left')
rdb_fre = tk.Radiobutton(frm1, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack(side='left')
var_boundary_cond.set(1)

# Cells
frm2 = ttk.Labelframe(root, relief="ridge", text="Cells parameters", labelanchor="n")
frm2.pack(side='left')
lbl_mass = tk.Label(frm2, text="mass")
lbl_mass.pack(side='left')
var_mass = tk.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_mass = tk.Spinbox(
    frm2, textvariable=var_mass, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_mass(var_mass.get()), width=6
    )
spn_mass.pack(side='left')
lbl_k = tk.Label(frm2, text=", k(Constant of springs between cells)")
lbl_k.pack(side='left')
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    frm2, textvariable=var_k, format="%.3f", from_=0.1, to=100., increment=0.001,
    command=lambda: change_k(var_k.get()), width=6
    )
spn_k.pack(side='left')

# Radio button
frm3 = ttk.Labelframe(root, relief="ridge", text="Click option", labelanchor="n")
frm3.pack(side='left')
var_radio_increase_decrease = tk.IntVar(root)
# Radio button 1st
r_increase = tk.Radiobutton(frm3, text="Check value", value=0, var=var_radio_increase_decrease)
r_increase.pack()
# Radio button 2st
r_increase = tk.Radiobutton(frm3, text="Increase", value=1, var=var_radio_increase_decrease)
r_increase.pack()
# Radio button 3nd
r_decrease = tk.Radiobutton(frm3, text="Decrease", value=2, var=var_radio_increase_decrease)
r_decrease.pack()
var_radio_increase_decrease.set(0)  # set default
# Label inc/dec step
lbl_step = tk.Label(frm3, text="Inc/dec step")
lbl_step.pack()
var_step = tk.StringVar(root)  # variable for spinbox-value
var_step.set(step_inc_dec)  # Initial value
spn_step = tk.Spinbox(
    frm3, textvariable=var_step, format="%.1f", from_=1, to=10, increment=1,
    command=lambda: change_inc_dec_step(var_step.get()), width=6
    )
spn_step.pack(side='left')

lbl_value_cell = tk.Label(root, text=" Cell(x0, y0)=0")
lbl_value_cell.pack(side='left')

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=100)
root.bind('<Configure>', on_change_window)
root.mainloop()
