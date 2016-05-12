#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
from math import exp, e
from tkinter import Tk, Frame, Label, CENTER, Scale, GROOVE, HORIZONTAL, LEFT, Button
import numpy as np
import matplotlib.pyplot as plt

from methods import runge_kutta, explicit_euler, recount_euler, get_next_m, progon

N = 8
L = 2 + 0.1 * N


def f(x, y, d_y):
    return y + L*x*(1 - x) + 2 * L + 2


def real_func(x):
    return -2 + exp(-x) + exp(x) + L*x*(-1 + x)


def main():
    from bokeh.plotting import figure, output_file, show
    # Anna Mokrushina & RomanDubinin drawer

    n = 30  # int(sys.argv[1])
    x_points = np.linspace(0, 1, n+1, endpoint=True)
    y_0 = 0
    right_boundary = e - 1.0/e + L
    h = 1.0/n
    colors = ["blue", "purple", "orange", "green"]

    methods = [
        {"f": explicit_euler, "name": "euler"},
        {"f": recount_euler,  "name": "euler with recount"},
        {"f": runge_kutta,    "name": "Runge-Kutta"}
    ]
    p = figure(title="", plot_width=1100, plot_height=600)
    y_points = list(map(real_func, x_points))
    p.line(x=x_points, y=y_points, color="black", legend="real", line_width=2)

    y_points = progon(x_points, right_boundary, h, L)
    p.line(x=x_points, y=y_points, color="red", legend="progon", line_width=2)

    for i in range(len(methods)):
        # y_points = real_func_method(x_points)

        last_m = 0
        y_points = [methods[i]["f"](f, x_points, y_0, last_m, h)]
        next_m = get_next_m(x_points, f, methods[i]["f"], y_0, right_boundary, last_m, h)

        while abs(next_m - last_m) > 0.01:
            y_points.append(methods[i]["f"](f, x_points, y_0, next_m, h))
            last_m = next_m
            next_m = get_next_m(x_points, f, methods[i]["f"], y_0, right_boundary, last_m, h)

        p.line(x=x_points, y=y_points[-1], color=colors[i], legend=methods[i]["name"])

    output_file("{0}_points.html".format(n))
    p.legend.orientation = "top_left"
    show(p)


def main2():
    methods = [
        (explicit_euler, 'r8-'),
        (recount_euler, 'rs-'),
        (runge_kutta, 'yd-'),
    ]

    def draw_all(n):
        plt.close()
        x_points = np.linspace(0, 1, n, endpoint=True)
        y_0 = 0
        right_boundary = e - 1.0/e + L

        for method, color in methods:
            last_m = 0
            y_points = method(x_points, 1.0 / n, y_0, last_m, f)
            next_m = get_next_m(x_points, f, method, y_0, right_boundary, last_m, 1.0 / n)

            while abs(next_m - last_m) > 0.01:
                y_points = method(x_points, 1.0 / n, y_0, next_m, f)
                last_m = next_m
                next_m = get_next_m(x_points, f, method, y_0, right_boundary, last_m, 1.0 / n)

            plt.plot(x_points, y_points, color, label=method.__doc__)
        plt.plot(x_points, list(map(real_func, x_points)), 'b*-', label='Original')
        plt.plot(x_points, progon(x_points, right_boundary, 1.0 / n, L), "r*-", label="progon")

        legend = plt.legend(loc='upper left', shadow=True, fontsize='x-large')
        legend.get_frame().set_facecolor('#00FFCC')

        mng = plt.get_current_fig_manager()

        mng.window.state('zoomed')
        # mng.full_screen_toggle()
        plt.show()

    window = Tk()
    # window.iconbitmap('Nik.ico')
    window.title('Селезнев Никита, ФИИТ-401')
    window.resizable(width=False, height=False)
    color1 = 'PaleGoldenrod'
    color2 = 'lightyellow'

    frame = Frame(window, bg=color1)
    frame.pack()

    l1 = Label(frame, text="Лабораторная работа №2.\n14 вариант.",
               justify=CENTER, font=("Helvetica", 12), bd=10, bg=color1)
    l1.pack()

    l = Label(frame, text="""f(xcdcdвсвыамывмывмывscsdcsd^2 * (x - 0.4),\ty0 = 0.5
    Решение численными методами
        1) Эйлера (явный)
        2) Эйлера (с пересчётом)
        3) Рунге-Кутта

        ?) Прогонка?
    а также точное решение.""", justify=LEFT, font=("Helvetica", 12), bd=10, bg=color2)
    l.pack()

    l2 = Label(frame, text="\nВыберите количество точек:", justify=CENTER, font=("Helvetica", 12), bd=0, bg=color1)
    l2.pack()

    w = Scale(frame, from_=20, to=500, resolution=10, length=300, bg=color1, borderwidth=0,
              relief=GROOVE, orient=HORIZONTAL, highlightthickness=0)
    w.pack()

    Label(frame, text='\n', bg=color1).pack()

    button = Button(frame, text="Нарисовать график",
                    font=("Helvetica", 12),
                    bg=color2,
                    command=lambda: draw_all(int(w.get())))
    button.pack()

    window.mainloop()


if __name__ == '__main__':
    main2()
