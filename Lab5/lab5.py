#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
from math import exp, e
from tkinter import Tk, Frame, Label, CENTER, Scale, GROOVE, HORIZONTAL, LEFT, Button
import numpy as np
import matplotlib.pyplot as plt

from methods import runge_kutta, explicit_euler, recount_euler, get_next_m, progon

N = 13
L = 2 + 0.1 * N


def f(x, y, d_y):
    return y + L*x*(1 - x) + 2*L + 2


def real_func(x):
    return -2 + exp(-x) + exp(x) + L*x*(-1 + x)


def draw_all(n):
        plt.close()
        x_points = np.linspace(0, 1, n, endpoint=True)
        y_0 = 0.0
        right_boundary = e - 1.0/e + L

        for method, color in [
            (explicit_euler, 'r8-'),
            (recount_euler, 'rs-'),
            (runge_kutta, 'yd-'),
        ]:
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


def main():
    window = Tk()
    # window.iconbitmap('Nik.ico')
    window.title('Селезнев Никита, ФИИТ-401')
    window.resizable(width=False, height=False)
    color1 = 'PaleGoldenrod'
    color2 = 'lightyellow'

    frame = Frame(window, bg=color1)
    frame.pack()

    Label(frame, text="Лабораторная работа №2.\n14 вариант.",
          justify=CENTER, font=("Helvetica", 12), bd=10, bg=color1).pack()

    Label(frame, text="""f(x) = y + Lx(1-x) + 2L + 2, L = 3.3,\ty0 = 0
    Решение численными методами
        1) Эйлера (явный)
        2) Эйлера (с пересчётом)
        3) Рунге-Кутта
        4) Прогонка
    а также точное решение.""", justify=LEFT, font=("Helvetica", 12), bd=10, bg=color2).pack()

    Label(frame, text="\nВыберите количество точек:", justify=CENTER, font=("Helvetica", 12), bd=0, bg=color1).pack()

    w = Scale(frame, from_=5, to=50, resolution=1, length=300, bg=color1, borderwidth=0,
              relief=GROOVE, orient=HORIZONTAL, highlightthickness=0)
    w.pack()

    Label(frame, text='\n', bg=color1).pack()

    Button(frame, text="Нарисовать график",
           font=("Helvetica", 12),
           bg=color2,
           command=lambda: draw_all(int(w.get()))).pack()

    window.mainloop()


if __name__ == '__main__':
    main()
