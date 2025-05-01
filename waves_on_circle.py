""" Waves on circle """
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d

""" Global variables """

""" Animation control """
is_play = False

""" Other parameters """
sigma_gauss = 0.03
amplitude_gauss = 0.1
phase_gauss_deg = 0.

wn_cos = 40.
amplitude_cos = 0.4
phase_cos_deg = 0.

scale_displacement = 1.
scale_velocity = 1.

""" Create figure and axes """
title_ax0 = "Waves on circle"
title_tk = title_ax0

x_min = -2.
x_max = 2.
y_min = -2.
y_max = 2.
z_min = -2.
z_max = 2.

fig = Figure()
ax0 = fig.add_subplot(121, projection='3d')
ax0.set_box_aspect((1, 1, 1))
ax0.grid()
ax0.set_title(title_ax0)
ax0.set_xlabel("x")
ax0.set_ylabel("y (velocity)")
ax0.set_zlabel("z")
ax0.set_xlim(x_min, x_max)
ax0.set_ylim(y_min, y_max)
ax0.set_zlim(z_min, z_max)

ax1 = fig.add_subplot(122)
ax1.set_title("Displacement and velocity")
ax1.set_xlabel("x")
ax1.set_ylabel("z")
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(z_min, z_max)
ax1.set_aspect("equal")
# ax1.set_aspect(30)
ax1.grid()


""" Embed in Tkinter """
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill="both")

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

""" Global objects of Tkinter """
var_phase_stp = tk.StringVar(root)
var_type = tk.IntVar(root)

var_sigma_gauss = tk.StringVar(root)
var_amplitude_gauss = tk.StringVar(root)
var_phase_gauss = tk.StringVar(root)

var_wn_cos = tk.StringVar(root)
var_amplitude_cos = tk.StringVar(root)
var_phase_cos = tk.StringVar(root)

""" Classes and functions """


class Counter:
    def __init__(self, is3d=None, ax=None, xy=None, z=None, label=""):
        self.is3d = is3d if is3d is not None else False
        self.ax = ax
        self.x, self.y = xy[0], xy[1]
        self.z = z if z is not None else 0
        self.label = label

        self.count = 0

        if not is3d:
            self.txt_step = self.ax.text(self.x, self.y, self.label + str(self.count))
        else:
            self.txt_step = self.ax.text2D(self.x, self.y, self.label + str(self.count))
            self.xz, self.yz, _ = proj3d.proj_transform(self.x, self.y, self.z, self.ax.get_proj())
            self.txt_step.set_position((self.xz, self.yz))

    def count_up(self):
        self.count += 1
        self.txt_step.set_text(self.label + str(self.count))

    def reset(self):
        self.count = 0
        self.txt_step.set_text(self.label + str(self.count))

    def get(self):
        return self.count


class WaveOnCircle:
    def __init__(self, ax, ax3d):
        self.ax = ax
        self.ax3d = ax3d

        self.theta = np.arange(0, 2 * np.pi, 0.01)
        self.displacement = self.theta * 0.
        self.displacement_buffer = self.theta * 0.
        self.velocity = self.theta * 0.

        self.r = 1

        self.x_d = self.r * np.cos(self.theta)
        self.y_d = self.r * np.sin(self.theta)

        self.circle_displacement, = self.ax.plot(self.x_d, self.y_d, color="red", linestyle="-", linewidth=2,
                                                 label="Displacement", alpha=0.5)

        self.x_v = self.r * np.cos(self.theta)
        self.y_v = self.r * np.sin(self.theta)
        self.circle_velocity, = self.ax.plot(self.x_v, self.y_v, color="green", linestyle="--", linewidth=1,
                                             label="Velocity", alpha=1)

        self.type = 1

        self.sigma_gauss = 0.05
        self.phase_gauss_deg = 0.
        self.amplitude_gauss = 0.05

        self.wn_cos = 20.
        self.phase_cos_deg = 0.
        self.amplitude_cos = 0.06

        self.k = 1.
        self.mass = 1.5

        self.scale_displacement = 1.
        self.scale_velocity = 1.

        self.x3d = self.x_d
        self.y3d = self.theta * 0.
        self.z3d = self.y_d
        self.circle_3d, = self.ax3d.plot(np.array(self.x3d), np.array(self.y3d), np.array(self.z3d),
                                         color="red", linewidth=1)

    def set_type(self, value):
        self.type = value
        if self.type == 1:
            self.set_gaussian()
        elif self.type == 2:
            self.set_cos()
        else:
            pass

    def _set_plot(self):
        self.x_d = (self.r + self.scale_displacement * self.displacement) * np.cos(self.theta)
        self.y_d = (self.r + self.scale_displacement * self.displacement) * np.sin(self.theta)

        self.circle_displacement.set_data(self.x_d, self.y_d)

        self.x_v = (self.r + self.scale_velocity * self.velocity) * np.cos(self.theta)
        self.y_v = (self.r + self.scale_velocity * self.velocity) * np.sin(self.theta)

        self.circle_velocity.set_data(self.x_v, self.y_v)

        self.circle_3d.set_xdata(np.array(self.x_d))
        self.velocity3d = self.scale_velocity * self.velocity
        self.circle_3d.set_ydata(np.array(self.velocity3d))
        self.circle_3d.set_3d_properties(np.array(self.y_d))

    def set_gaussian(self):
        displacement_center = (self.amplitude_gauss * 1 / np.sqrt(2 * np.pi * self.sigma_gauss ** 2) *
                               np.e ** (- (self.theta - np.deg2rad(180)) ** 2 / (2 * self.sigma_gauss ** 2)))
        index_roll = int(((len(self.theta) - 0) / (2 * np.pi)) * (np.deg2rad(self.phase_gauss_deg) - np.pi))
        self.displacement = np.roll(displacement_center, index_roll)

        self._set_plot()

    def set_sigma_gauss(self, value):
        self.sigma_gauss = value
        if self.type == 1:
            self.set_gaussian()

    def set_phase_gauss_deg(self, value):
        self.phase_gauss_deg = value
        if self.type == 1:
            self.set_gaussian()

    def set_amplitude_gauss(self, value):
        self.amplitude_gauss = value
        if self.type == 1:
            self.set_gaussian()

    def set_cos(self):
        self.displacement = self.amplitude_cos * np.cos(self.wn_cos * (self.theta - np.deg2rad(self.phase_cos_deg)))
        self._set_plot()

    def set_wn_cos(self, value):
        self.wn_cos = value
        if self.type == 2:
            self.set_cos()

    def set_phase_cos_deg(self, value):
        self.phase_cos_deg = value
        if self.type == 2:
            self.set_cos()

    def set_amplitude_cos(self, value):
        self.amplitude_cos = value
        if self.type == 2:
            self.set_cos()

    def reset(self):
        self.displacement = self.theta * 0.
        self.displacement_buffer = self.theta * 0.
        self.velocity = self.theta * 0.

        self._set_plot()

    def update_wave(self):
        for i in range(len(self.theta)):
            if i == 0:
                dd = - (self.displacement[-1] - self.displacement[i]) - (self.displacement[i + 1] - self.displacement[i])
            elif i == len(self.theta) - 1:
                dd = - (self.displacement[i - 1] - self.displacement[i]) - (self.displacement[0] - self.displacement[i])
            else:
                dd = - (self.displacement[i - 1] - self.displacement[i]) - (self.displacement[i + 1] - self.displacement[i])
            force = - self.k * dd
            a = force / self.mass
            self.velocity[i] = self.velocity[i] + a
            self.displacement_buffer[i] = self.displacement_buffer[i] + self.velocity[i]

        self.displacement = self.displacement_buffer.copy()
        self._set_plot()

    def adjust_scale(self):
        max_d = np.max(abs(self.displacement))
        max_v = np.max(abs(self.velocity))

        if max_v != 0:
            self.scale_velocity = max_d / max_v


def create_parameter_setter():
    # Type of wave
    frm_type = ttk.Labelframe(root, relief="ridge", text="Type", labelanchor='n')
    frm_type.pack(side="left", fill=tk.Y)

    # var_type = tk.IntVar(root)
    rd_type_gauss = tk.Radiobutton(frm_type, text="Gaussian", value=1, variable=var_type,
                                   command=lambda: wave_on_circle.set_type(var_type.get()))
    rd_type_gauss.pack(anchor=tk.W)
    rd_type_cos = tk.Radiobutton(frm_type, text="Cosine", value=2, variable=var_type,
                                 command=lambda: wave_on_circle.set_type(var_type.get()))
    rd_type_cos.pack(anchor=tk.W)
    var_type.set(1)

    # Parameter of gaussian
    frm_gauss = ttk.Labelframe(root, relief="ridge", text="Gaussian (initial)", labelanchor='n')
    frm_gauss.pack(side="left", fill=tk.Y)

    lbl_sigma = tk.Label(frm_gauss, text="Sigma")
    lbl_sigma.pack(side='left')
    # var_sigma_gauss = tk.StringVar(root)
    var_sigma_gauss.set(str(sigma_gauss))
    spn_sigma = tk.Spinbox(
        frm_gauss, textvariable=var_sigma_gauss, format="%.2f", from_=0.01, to=1.0, increment=0.01,
        command=lambda: wave_on_circle.set_sigma_gauss(float(var_sigma_gauss.get())), width=5
    )
    spn_sigma.pack(side="left")

    lbl_phase_gauss = tk.Label(frm_gauss, text="Phase (degree)")
    lbl_phase_gauss.pack(side='left')
    # var_phase_gauss = tk.StringVar(root)
    var_phase_gauss.set(str(phase_gauss_deg))
    spn_phase_gauss = tk.Spinbox(
        frm_gauss, textvariable=var_phase_gauss, format="%.1f", from_=-360, to=360, increment=1,
        command=lambda: wave_on_circle.set_phase_gauss_deg(float(var_phase_gauss.get())), width=5
    )
    spn_phase_gauss.pack(side="left")

    lbl_amplitude_gauss = tk.Label(frm_gauss, text="Amp.")
    lbl_amplitude_gauss.pack(side='left')
    # var_amplitude_gauss = tk.StringVar(root)
    var_amplitude_gauss.set(str(amplitude_gauss))
    spn_amplitude_gauss = tk.Spinbox(
        frm_gauss, textvariable=var_amplitude_gauss, format="%.2f", from_=0.01, to=5.0, increment=0.01,
        command=lambda: wave_on_circle.set_amplitude_gauss(float(var_amplitude_gauss.get())), width=5
    )
    spn_amplitude_gauss.pack(side="left")

    # Parameter of cosine
    frm_cos = ttk.Labelframe(root, relief="ridge", text="Cosine (initial)", labelanchor='n')
    frm_cos.pack(side="left", fill=tk.Y)

    lbl_wn = tk.Label(frm_cos, text="Wave number")
    lbl_wn.pack(side='left')
    # var_wn_coss = tk.StringVar(root)
    var_wn_cos.set(str(wn_cos))
    spn_wn_cos = tk.Spinbox(
        frm_cos, textvariable=var_wn_cos, format="%.1f", from_=1, to=50, increment=1,
        command=lambda: wave_on_circle.set_wn_cos(float(var_wn_cos.get())), width=5
    )
    spn_wn_cos.pack(side="left")

    lbl_phase_cos = tk.Label(frm_cos, text="Phase (degree)")
    lbl_phase_cos.pack(side='left')
    # var_phase_cos = tk.StringVar(root)
    var_phase_cos.set(str(phase_cos_deg))
    spn_phase_cos = tk.Spinbox(
        frm_cos, textvariable=var_phase_cos, format="%.1f", from_=-360, to=360, increment=1,
        command=lambda: wave_on_circle.set_phase_cos_deg(float(var_phase_cos.get())), width=5
    )
    spn_phase_cos.pack(side="left")

    lbl_amplitude_cos = tk.Label(frm_cos, text="Amp.")
    lbl_amplitude_cos.pack(side='left')
    # var_amplitude_cos = tk.StringVar(root)
    var_amplitude_cos.set(str(amplitude_cos))
    spn_amplitude_cos = tk.Spinbox(
        frm_cos, textvariable=var_amplitude_cos, format="%.2f", from_=0.01, to=5.0, increment=0.01,
        command=lambda: wave_on_circle.set_amplitude_cos(float(var_amplitude_cos.get())), width=5
    )
    spn_amplitude_cos.pack(side="left")

    # Adjust scale
    frm_adjust = ttk.Labelframe(root, relief="ridge", text="Adjust scale", labelanchor='n')
    frm_adjust.pack(side="left", fill=tk.Y)

    btn_set_adjust = tk.Button(frm_adjust, text="Adjust", command=lambda: wave_on_circle.adjust_scale())
    btn_set_adjust.pack(fill=tk.X)


def create_animation_control():
    frm_anim = ttk.Labelframe(root, relief="ridge", text="Animation", labelanchor="n")
    frm_anim.pack(side="left", fill=tk.Y)
    btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
    btn_play.pack(fill=tk.X)
    btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
    btn_reset.pack(fill=tk.X)
    # btn_clear = tk.Button(frm_anim, text="Clear path", command=lambda: aaa))
    # btn_clear.pack(fill=tk.X)


def create_center_lines():
    line_axis_x = art3d.Line3D([0., 0.], [0., 0.], [z_min, z_max], color="gray", ls="-.", linewidth=1)
    ax0.add_line(line_axis_x)
    line_axis_y = art3d.Line3D([x_min, x_max], [0., 0.], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(line_axis_y)
    line_axis_z = art3d.Line3D([0., 0.], [y_min, y_max], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(line_axis_z)


def create_circle(ax, x, y, z, z_dir, edge_col, fill_flag, line_width, line_style, label):
    if label != "":
        c_spin_axis_guide = Circle((x, y), 1., ec=edge_col, fill=fill_flag,
                                   linewidth=line_width, linestyle=line_style, label=label)
    else:
        c_spin_axis_guide = Circle((x, y), 1., ec=edge_col, fill=fill_flag,
                                   linewidth=line_width, linestyle=line_style)
    ax.add_patch(c_spin_axis_guide)
    art3d.pathpatch_2d_to_3d(c_spin_axis_guide, z=z, zdir=z_dir)


def draw_static_diagrams():
    create_center_lines()
    # create_circle(ax0, 0., 0., 0., "x", "gray", False, 0.5, "--", "")
    create_circle(ax0, 0., 0., 0., "y", "gray", False, 0.5,"--", "")
    # create_circle(ax0, 0., 0., 0., "z", "gray", False, 0.5, "--", "")


def reset():
    global is_play
    is_play = False
    cnt.reset()
    wave_on_circle.reset()


def switch():
    global is_play
    is_play = not is_play


def update(f):
    if is_play:
        cnt.count_up()
        wave_on_circle.update_wave()


""" main loop """
if __name__ == "__main__":
    cnt = Counter(ax=ax0, is3d=True, xy=np.array([x_min, y_max]), z=z_max, label="Step=")
    draw_static_diagrams()
    create_animation_control()
    create_parameter_setter()

    wave_on_circle = WaveOnCircle(ax1, ax0)

    wave_on_circle.set_sigma_gauss(sigma_gauss)
    wave_on_circle.set_amplitude_gauss(amplitude_gauss)
    wave_on_circle.set_phase_gauss_deg(phase_gauss_deg)

    wave_on_circle.set_wn_cos(wn_cos)
    wave_on_circle.set_amplitude_cos(amplitude_cos)
    wave_on_circle.set_phase_cos_deg(phase_cos_deg)

    ax1.legend(loc='lower right', fontsize=8)

    anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
    root.mainloop()
