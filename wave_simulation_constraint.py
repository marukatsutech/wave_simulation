# Wave simulation with constraint
import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.animation as animation
import numpy as np
import matplotlib.patches as patches


def change_constraint_mode():
    global mode_constraint
    mode_constraint = var_mode_constraint.get()


def change_ratio_reduction_constraint(value):
    global ratio_reduction
    ratio_reduction = float(value)


def change_width_constraint(value):
    global width_constraint
    width_constraint = float(value)


def change_k(value):
    global k
    k = float(value)


def change_mass(value):
    global mass
    mass = float(value)


def reset_string():
    global y, v
    y = y.copy() * 0.
    v = v.copy() * 0.
    string1.set_ydata(y)


def change_sigma(value):
    global sigma
    sigma = float(value)


def impact_gauss():
    global y
    gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- (x - impact_x) ** 2 / (2 * sigma ** 2))
    y = y.copy() + gaussian
    if var_boundary_cond.get() == 1:
        y[0] = 0.
        y[len(x) - 1] = 0.
    string1.set_ydata(y)


def update_constraint():
    global x_constraint1, x_constraint2, tx_area1, tx_area2
    # Areas
    area_constraint1.set_width(width_constraint)
    area_constraint1.set_x(x_constraint1 - width_constraint / 2.)
    area_constraint2.set_width(width_constraint)
    area_constraint2.set_x(x_constraint2 - width_constraint / 2.)
    # lines
    x_constraint1 = x_min + float(scale_var_x_constraint1.get()) * (x_max - x_min) / 100.
    line_constraint1.set_data([x_constraint1, x_constraint1], [y_min, y_max])
    x_constraint2 = x_min + float(scale_var_x_constraint2.get()) * (x_max - x_min) / 100.
    line_constraint2.set_data([x_constraint2, x_constraint2], [y_min, y_max])
    # Texts
    tx_area1.set_position((x_constraint1, y_max * 0.9))
    tx_area2.set_position((x_constraint2, y_max * 0.9))


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    global cnt, tx_step, line_impact, impact_x, string1, y, y_buffer, v
    # Update items
    # Step
    tx_step.set_text(" Step(as t)=" + str(cnt))
    # Impact line
    impact_x = x_min + float(scale_var.get()) * (x_max - x_min) / 100.
    line_impact.set_data([impact_x, impact_x], [y_min, y_max])
    # Constraint
    update_constraint()
    if is_play:
        # String
        for i in range(len(x)):
            if var_boundary_cond.get() == 1:
                if i == 0:
                    dy = 0.
                elif i == 1:
                    dy = - (0. - y[i]) - (y[i + 1] - y[i])
                elif i == len(x) - 2:
                    dy = - (y[i - 1] - y[i]) - (0. - y[i])
                elif i == len(x) - 1:
                    dy = 0.
                else:
                    dy = - (y[i - 1] - y[i]) - (y[i + 1] - y[i])
            else:
                if i == 0:
                    dy = - (y[i + 1] - y[i])
                elif i == len(x) - 1:
                    dy = - (y[i - 1] - y[i])
                else:
                    dy = - (y[i - 1] - y[i]) - (y[i + 1] - y[i])
            x_i = i * (x_max - x_min) / (num_of_mass - 1) + x_min
            if np.abs(x_i - x_constraint1) < width_constraint / 2. or np.abs(x_i - x_constraint2) < width_constraint / 2.:
                force = - k * dy
                a = force / mass
                if mode_constraint == 0:
                    v[i] = v[i] + a * ratio_reduction
                    y_buffer[i] = y[i] + v[i]
                elif mode_constraint == 1:
                    v[i] = (v[i] + a) * ratio_reduction
                    y_buffer[i] = y[i] + v[i]
                else:
                    v[i] = v[i] + a
                    y_buffer[i] = (y[i] + v[i]) * ratio_reduction
            else:
                force = - k * dy
                a = force / mass
                v[i] = v[i] + a
                y_buffer[i] = y[i] + v[i]
        y = y_buffer.copy()
        string1.set_ydata(y)
        cnt += 1


# Global variables
x_min = 0.
x_max = 20.
y_min = -4.
y_max = 4.

num_of_mass = 400
mass = 10.
mass_init = mass
k = 10.
k_init = k

impact_x = 0.
sigma = 0.2
sigma_init = sigma

k_constraint = 10.
x_constraint1 = 5.
x_constraint2 = 15.
width_constraint = 1.0
height_constraint = y_max - y_min

ratio_reduction = 1.
mode_constraint = 0

# For UI
is_play = False
cnt = 0

# Generate tkinter
root = tkinter.Tk()
root.title("Wave simulation with constraint")

# Generate figure and axes
fig = Figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_title('Wave simulation with constraint')
ax1.set_ylabel('y')
ax1.grid()

# Generate items Note;  Some objects of ax.plot need ',' to use "set_* method".
tx_step = ax1.text(x_min, y_max * 0.9, ' Step=' + str(0) + ',k=' + str(0) + ',omega=' + str(0) + '/step')
line_impact, = ax1.plot([impact_x, impact_x], [y_min, y_max], linestyle='-.', c='red')
x = np.linspace(x_min, x_max, num_of_mass)
y = x * 0.
y_buffer = y
v = x * 0.
string1, = ax1.plot(x, y, linestyle='-')

area_constraint1 = patches.Rectangle(
    xy=(x_constraint1 - width_constraint / 2., y_min),
    width=width_constraint, height=height_constraint, fc='blue', alpha=0.2
)
ax1.add_patch(area_constraint1)
area_constraint2 = patches.Rectangle(
    xy=(x_constraint2 - width_constraint / 2., y_min),
    width=width_constraint, height=height_constraint, fc='blue', alpha=0.2
)
ax1.add_patch(area_constraint2)
line_constraint1, = ax1.plot([x_constraint1, x_constraint1], [y_min, y_max], linestyle='-.', c='blue')
line_constraint2, = ax1.plot([x_constraint2, x_constraint2], [y_min, y_max], linestyle='-.', c='blue')
tx_area1 = ax1.text(x_constraint1, y_max * 0.9, 'Area1')
tx_area2 = ax1.text(x_constraint2, y_max * 0.9, 'Area2')

# Embed a figure in canvas
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack()

# Tkinter widgets
# Toolbar
toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

# Widgets to set parameters
# Boundary condition
frm_boundary = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
frm_boundary.pack(side='left')
var_boundary_cond = tkinter.IntVar(root)
rdb_fxe = tkinter.Radiobutton(frm_boundary, text="Fixed end, ", value=1, var=var_boundary_cond)
rdb_fxe.pack(side='left')
rdb_fre = tkinter.Radiobutton(frm_boundary, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack(side='left')
var_boundary_cond.set(1)

# String
frm_string = ttk.Labelframe(root, relief="ridge", text="String", labelanchor="n")
frm_string.pack(side='left')
lbl_mass = tkinter.Label(frm_string, text="mass")
lbl_mass.pack(side='left')
var_mass = tkinter.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_mass = tkinter.Spinbox(
    frm_string, textvariable=var_mass, format="%.1f", from_=1., to=20., increment=1.,
    command=lambda: change_mass(var_mass.get()), width=4
)
spn_mass.pack(side='left')
lbl_k = tkinter.Label(frm_string, text=", k")
lbl_k.pack(side='left')
var_k = tkinter.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tkinter.Spinbox(
    frm_string, textvariable=var_k, format="%.1f", from_=1., to=10., increment=1.,
    command=lambda: change_k(var_k.get()), width=4
)
spn_k.pack(side='left')

# Constraint
frm_constraint = ttk.Labelframe(root, relief="ridge", text="Constraint", labelanchor="n")
frm_constraint .pack(side='left')
var_mode_constraint = tkinter.IntVar(value=mode_constraint)
rdb0_mode_constraint = tkinter.Radiobutton(frm_constraint, text="k(=acceleration)", command=change_constraint_mode,
                                           variable=var_mode_constraint, value=0)
rdb1_mode_constraint = tkinter.Radiobutton(frm_constraint, text="velocity", command=change_constraint_mode,
                                           variable=var_mode_constraint, value=1)
rdb2_mode_constraint = tkinter.Radiobutton(frm_constraint, text="y", command=change_constraint_mode,
                                           variable=var_mode_constraint, value=2)
rdb0_mode_constraint.grid(column=0, row=0)
rdb1_mode_constraint.grid(column=1, row=0)
rdb2_mode_constraint.grid(column=0, row=1, sticky=tkinter.W)
lbl_reduction_ratio_constraint = tkinter.Label(frm_constraint, text="Reduction ratio")
lbl_reduction_ratio_constraint.grid(column=0, row=2)
var_reduction_ratio_constraint = tkinter.StringVar(root)  # variable for spinbox-value
var_reduction_ratio_constraint.set(ratio_reduction)  # Initial value
spn_reduction_ratio_constraint = tkinter.Spinbox(
    frm_constraint, textvariable=var_reduction_ratio_constraint, format="%.2f", from_=0.0, to=1., increment=0.01,
    command=lambda: change_ratio_reduction_constraint(var_reduction_ratio_constraint.get()), width=4
)
spn_reduction_ratio_constraint.grid(column=1, row=2, sticky=tkinter.W)
lbl_width_constraint = tkinter.Label(frm_constraint, text="Width")
lbl_width_constraint.grid(column=0, row=3)
var_width_constraint = tkinter.StringVar(root)  # variable for spinbox-value
var_width_constraint.set(width_constraint)  # Initial value
spn_width_constraint = tkinter.Spinbox(
    frm_constraint, textvariable=var_width_constraint, format="%.2f", from_=0.1, to=5., increment=0.1,
    command=lambda: change_width_constraint(var_width_constraint.get()), width=4
)
spn_width_constraint.grid(column=1, row=3, sticky=tkinter.W)
lbl_x_constraint1 = tkinter.Label(frm_constraint, text="Area1")
lbl_x_constraint1.grid(column=0, row=4)
scale_var_x_constraint1 = tkinter.StringVar(root)
s_x_constraint1 = tkinter.Scale(frm_constraint, variable=scale_var_x_constraint1, orient='horizontal', length=200)
s_x_constraint1.grid(column=1, row=4)
scale_var_x_constraint1.set(25)
lbl_x_constraint2 = tkinter.Label(frm_constraint, text="Area2")
lbl_x_constraint2.grid(column=0, row=5)
scale_var_x_constraint2 = tkinter.StringVar(root)
s_x_constraint2 = tkinter.Scale(frm_constraint, variable=scale_var_x_constraint2, orient='horizontal', length=200)
s_x_constraint2.grid(column=1, row=5)
scale_var_x_constraint2.set(75)

# Gaussian curve
frm_gaussian = ttk.Labelframe(root, relief="ridge", text="Gaussian curve", labelanchor="n")
frm_gaussian .pack(side='left')
lbl_sigma = tkinter.Label(frm_gaussian, text="Sigma")
lbl_sigma.pack(side='left')
var_sgm = tkinter.StringVar(root)  # variable for spinbox-value
var_sgm.set(sigma_init)  # Initial value
spn_wn = tkinter.Spinbox(
    frm_gaussian, textvariable=var_sgm, format="%.1f", from_=0.1, to=1.0, increment=0.1,
    command=lambda: change_sigma(var_sgm.get()), width=4
)
spn_wn.pack(side='left')
lbl_position = tkinter.Label(frm_gaussian, text="Position")
lbl_position.pack(side='left')
scale_var = tkinter.StringVar(root)
s = tkinter.Scale(frm_gaussian, variable=scale_var, orient='horizontal', length=200)
s.pack(side='left')
scale_var.set(50)
btn_gauss = tkinter.Button(frm_gaussian, text="Push", command=impact_gauss)
btn_gauss .pack(side='left')

# Play/pause
btn_play_pause = tkinter.Button(root, text="Play/Pause", command=switch)
btn_play_pause.pack(side='left')

# Reset (flat)
lbl_space = tkinter.Label(root, text=" ")
lbl_space.pack(side='left')
btn_reset = tkinter.Button(root, text="Flat", command=reset_string)
btn_reset .pack(side='left')

# main loop
anim = animation.FuncAnimation(fig, update, interval=50)
tkinter.mainloop()
