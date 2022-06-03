# Wave simulation 1D
import tkinter
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.animation as animation
import numpy as np


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


def change_sigma(value):
    global sigma
    sigma = float(value)


def change_wn(value):
    global wave_number
    wave_number = float(value)


def impact_sine():
    global y
    if var_boundary_cond.get() == 1:
        sine = np.sin(wave_number * x * np.pi / (x_max - x_min))
        y = y.copy() + sine
        y[0] = 0.
        y[len(x) - 1] = 0.
    else:
        sine = np.sin(wave_number * x * np.pi / (x_max - x_min) - (x_max - x_min))
        y = y.copy() + sine
        print('test')


def impact_gauss():
    global y
    gaussian = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.e ** (- (x - impact_x) ** 2 / (2 * sigma ** 2))
    y = y.copy() + gaussian
    if var_boundary_cond.get() == 1:
        y[0] = 0.
        y[len(x) - 1] = 0.


def update(f):
    global tx_step, line_impact, impact_x, string1, y, y_buffer, v
    # Update items
    # Step
    tx_step.set_text(" Step(as t)=" + str(f))
    # Impact line
    impact_x = x_min + float(scale_var.get()) * (x_max - x_min) / 100.
    line_impact.set_data([impact_x, impact_x], [y_min, y_max])
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
        force = - k * dy
        a = force / mass
        v[i] = v[i] + a
        y_buffer[i] = y[i] + v[i]
    y = y_buffer.copy()
    string1.set_ydata(y)


# Global variables
x_min = 0.
x_max = 8.
y_min = -4.
y_max = 4.

num_of_mass = 400
mass = 10.
mass_init = mass
k = 10.
k_init = k

impact_x = 0.
sigma = 0.1
sigma_init = sigma
wave_number = 1
wn_init = wave_number

# Generate tkinter
root = tkinter.Tk()
root.title("Wave simulation")

# Generate figure and axes
fig = Figure(figsize=(8, 4))
ax1 = fig.add_subplot(111)
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_title('Wave simulation')
ax1.set_ylabel('y')
ax1.grid()

# Generate items Note;  Some objects of ax.plot need ',' to use "set_* method".
tx_step = ax1.text(x_min, y_max * 0.9, ' Step=' + str(0) + ',k=' + str(0) + ',omega=' + str(0) + '/step')
line_impact, = ax1.plot([impact_x, impact_x], [y_min, y_max], linestyle='-.', c='red')
x = np.linspace(0, x_max, num_of_mass)
y = x * 0.
y_buffer = y
v = x * 0.
string1, = ax1.plot(x, y, linestyle='-', label="String")
tx_comment1 = ax1.text(x_min, y_min * 0.8, ' Number of point mass=' + str(num_of_mass))
tx_comment2 = ax1.text(x_min, y_min * 0.9, ' k;Spring constant between points')

# Embed a figure in canvas
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack()

# Animation
anim = animation.FuncAnimation(fig, update, interval=50)

# Tkinter widgets

# Toolbar
toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack(expand=True, fill='both')

# Widgets to set parameters

# Boundary condition
frm1 = ttk.Labelframe(root, relief="ridge", text="Boundary condition", labelanchor="n", width=100)
var_boundary_cond = tkinter.IntVar(root)
rdb_fxe = tkinter.Radiobutton(frm1, text="Fixed end, ", value=1, var=var_boundary_cond)
rdb_fxe.pack(side='left')
rdb_fre = tkinter.Radiobutton(frm1, text="Free end", value=2, var=var_boundary_cond)
rdb_fre.pack(side='left')
var_boundary_cond.set(1)
frm1.pack(side='left')

# String
frm2 = ttk.Labelframe(root, relief="ridge", text="String parameters", labelanchor="n")
lbl_mass = tkinter.Label(frm2, text="mass")
lbl_mass.pack(side='left')
var_mass = tkinter.StringVar(root)  # variable for spinbox-value
var_mass.set(mass_init)  # Initial value
spn_mass = tkinter.Spinbox(
    frm2, textvariable=var_mass, format="%.1f", from_=1., to=20., increment=1.,
    command=lambda: change_mass(var_mass.get()), width=4
    )
spn_mass.pack(side='left')
lbl_k = tkinter.Label(frm2, text=", k")
lbl_k.pack(side='left')
var_k = tkinter.StringVar(root)  # variable for spinbox-value
var_k.set(k_init)  # Initial value
spn_k = tkinter.Spinbox(
    frm2, textvariable=var_k, format="%.1f", from_=1., to=10., increment=1.,
    command=lambda: change_k(var_k.get()), width=4
    )
spn_k.pack(side='left')
frm2.pack(side='left')

# Gaussian curve
frm3 = ttk.Labelframe(root, relief="ridge", text="Gaussian curve", labelanchor="n")
lbl_sigma = tkinter.Label(frm3, text="sigma")
lbl_sigma.pack(side='left')
var_sgm = tkinter.StringVar(root)  # variable for spinbox-value
var_sgm.set(sigma_init)  # Initial value
spn_wn = tkinter.Spinbox(
    frm3, textvariable=var_sgm, format="%.1f", from_=0.1, to=1.0, increment=0.1,
    command=lambda: change_sigma(var_sgm.get()), width=4
    )
spn_wn.pack(side='left')
lbl_position = tkinter.Label(frm3, text=", position")
lbl_position.pack(side='left')
scale_var = tkinter.StringVar(root)
s = tkinter.Scale(frm3, variable=scale_var, orient='horizontal', length=100)
s.pack(side='left')
scale_var.set(50)
btn_gauss = tkinter.Button(frm3, text="Push", command=impact_gauss)
btn_gauss .pack(side='left')
frm3.pack(side='left')

# Sine curve
frm4 = ttk.Labelframe(root, relief="ridge", text="Sine curve", labelanchor="n")
lbl_sine = tkinter.Label(frm4, text="wave number")
lbl_sine.pack(side='left')
var_wn = tkinter.StringVar(root)  # variable for spinbox-value
var_wn.set(wn_init)  # Initial value
spn_wn = tkinter.Spinbox(
    frm4, textvariable=var_wn, format="%.1f", from_=1, to=32, increment=1,
    command=lambda: change_wn(var_wn.get()), width=4
    )
spn_wn.pack(side='left')
btn_sine = tkinter.Button(frm4, text="Push", command=impact_sine)
btn_sine.pack(side='left')
frm4.pack(side='left')

# Reset (flat)
lbl_space = tkinter.Label(root, text=" ")
lbl_space.pack(side='left')
btn_reset = tkinter.Button(root, text="Flat", command=reset_string)
btn_reset .pack(side='left')

# main loop
tkinter.mainloop()
