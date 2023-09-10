# Wave simulation (Detail)

import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
import matplotlib.ticker as ticker


def change_k(value):
    global k
    k = float(value)


def change_mass(value):
    global mass
    mass = float(value)


def clear_cells():
    global is_play, cnt
    is_play = False
    cnt = 0
    for i in range(num_of_mass):
        y[i] = 0.
        dl[i] = 0.
        dr[i] = 0.
        f[i] = 0.
        a[i] = 0.
        v[i] = 0.
    update_cells()
    update_tree()


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
        if i == 0:
            dl[i] = y[i] - y_boundary_left
        else:
            dl[i] = y[i] - y[i - 1]
        if i == num_of_mass - 1:
            dr[i] = y[i] - y_boundary_right
        else:
            dr[i] = y[i] - y[i + 1]
        f[i] = - k * (dl[i] + dr[i])
        a[i] = f[i] / mass


def update_tree():
    values_y = ['y']
    values_dl = ['d_left']
    values_dr = ['d_right']
    values_f = ['f']
    values_a = ['a']
    values_v = ['v']
    for i in range(num_of_mass):
        values_y.append(y[i])
        values_dl.append(dl[i])
        values_dr.append(dr[i])
        values_f.append(f[i])
        values_a.append(f[i])
        values_v.append(v[i])
    tree.item(0, values=values_y)
    tree.item(1, values=values_dl)
    tree.item(2, values=values_dr)
    tree.item(3, values=values_f)
    tree.item(4, values=values_a)
    tree.item(5, values=values_v)


def update_cells():
    global lines_mass, lines_spring_left, lines_spring_right, lines_velocity
    # draw cells
    if var_boundary_cond.get() == 1:    # Fixed end
        y_boundary_left = 0.
        y_boundary_right = 0.
    else:                               # Free end
        y_boundary_left = y[0]
        y_boundary_right = y[num_of_mass - 1]
    for j in range(num_of_mass):
        lines_mass[j].set_data([j - 0.5, j + 0.5], [y[j], y[j]])
        lines_velocity[j].set_data([j - 0.5, j + 0.5], [v[j], v[j]])
        if j == 0:
            lines_spring_left[j].set_data([j - 0.5, j - 0.5], [y[j], y_boundary_left])
        else:
            lines_spring_left[j].set_data([j - 0.5, j - 0.5], [y[j], y[j - 1]])
        if j == num_of_mass - 1:
            lines_spring_right[j].set_data([j + 0.5, j + 0.5], [y[j], y_boundary_right])
        else:
            lines_spring_right[j].set_data([j + 0.5, j + 0.5], [y[j], y[j + 1]])


def mouse_motion(event):
    if event.dblclick == 1:
        # print("double click")
        pass
    elif event.button == 1:
        # print("left click")
        if str(event.xdata) != "None" and str(event.ydata) != "None":
            # print(event.xdata, event.ydata)
            if 0 <= round(event.xdata) <= num_of_mass - 1:
                if var_radio_y_or_v.get() == 0:
                    y[round(event.xdata)] = round(event.ydata)
                else:
                    v[round(event.xdata)] = round(event.ydata)
                update_cells()
                eval_cells()
                update_tree()
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
    update_tree()


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
        update_tree()


# Global variables
x_min = - 2.
x_max = 34.
y_min = - 10.
y_max = 10.

# Animation control
cnt = 0
is_play = False

# Parameter
num_of_mass = 32
mass = 10.
mass_init = mass
k = 1.
k_init = k

x = np.linspace(0, x_max, num_of_mass)
y = x * 0.
dl = x * 0.
dr = x * 0.
f = x * 0.
a = x * 0.
v = x * 0.

# Generate figure and axes
fig = Figure()
ax = fig.add_subplot(111)
ax.set_title("Wave simulation (Detail)")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
x_ticks = []
for ix in range(num_of_mass):
    x_ticks.append(ix)
ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
y_ticks = []
for iy in range(int(y_max-y_min)):
    y_ticks.append(y_min + iy)
ax.yaxis.set_major_locator(ticker.FixedLocator(y_ticks))


# ax.set_aspect("equal")
ax.grid()
# ax.invert_yaxis()

# Generate graphic items
x_cells0 = []
y_cells0 = []
s_cells0 = []  # maker size

tx_step = ax.text(x_min, y_max * 0.95, "Step=" + str(0))

lines_mass = []
lines_spring_left = []
lines_spring_right = []
lines_velocity = []
for iii in range(num_of_mass):
    line_m, = ax.plot([iii - 0.5, iii + 0.5], [0., 0.], linewidth=6)
    lines_mass.append(line_m)
    line_sl, = ax.plot([iii - 0.5, iii - 0.5], [0., 0.], linewidth=1, color='gray')
    lines_spring_left.append(line_sl)
    line_sr, = ax.plot([iii + 0.5, iii + 0.5], [0., 0.], linewidth=1, color='gray')
    lines_spring_right.append(line_sr)
    line_v, = ax.plot([iii - 0.5, iii + 0.5], [0., 0.], linewidth=1, color='blue', linestyle='--')
    lines_velocity.append(line_v)


# Tkinter
root = tk.Tk()
root.title("Wave simulation (Digitized detail)")
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
    command=lambda: change_mass(var_mass.get()), width=5
    )
spn_mass.pack(side='left')
lbl_k = tk.Label(frm2, text=", k(Constant of springs between cells)")
lbl_k.pack(side='left')
var_k = tk.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tk.Spinbox(
    frm2, textvariable=var_k, format="%.1f", from_=1., to=100., increment=1.,
    command=lambda: change_k(var_k.get()), width=5
    )
spn_k.pack(side='left')

# Radio button
frm3 = ttk.Labelframe(root, relief="ridge", text="Click option", labelanchor="n")
frm3.pack(side='left')
var_radio_y_or_v = tk.IntVar(root)
# Radio button 1st
r_y = tk.Radiobutton(frm3, text="y", value=0, var=var_radio_y_or_v)
r_y.pack()
# Radio button 2nd
r_v = tk.Radiobutton(frm3, text="v", value=1, var=var_radio_y_or_v)
r_v.pack()
var_radio_y_or_v.set(0)  # set default


# Tree view
column = ['variable']
for mm in range(num_of_mass):
    column.append('mass'+str(mm))

tree = ttk.Treeview(root, columns=column, height=6)

tree.column('variable', anchor='center', width=32, stretch='no')
for ii in range(num_of_mass):
    tree.column('mass'+str(ii), anchor='center', width=32, stretch='no')

tree.heading('variable', text='variable', anchor='center')
for jj in range(num_of_mass):
    tree.heading('mass'+str(jj), text='m'+str(jj), anchor='center')

dummy_value = ['variable']
for nn in range(num_of_mass):
    dummy_value.append(nn)

for kk in range(6):
    tree.insert(parent='', index='end', iid=kk, values=dummy_value)

tree.pack(side='left')

update_cells()
update_tree()

# Draw animation
anim = animation.FuncAnimation(fig, update, interval=100)
root.bind('<Configure>', on_change_window)
root.mainloop()
