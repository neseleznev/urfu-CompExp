#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
from tkinter import Tk, Frame, Button, Scale, HORIZONTAL, Label, LEFT, CENTER, FLAT, GROOVE, filedialog

from methods import explicit_euler, euler_with_recount, cauchy, \
    runge_kutta, implicit_euler, taylor_3, taylor_4, extrapolation_adams

GROUP = 3

tasks = {
    '1': {
        'y0': 0.1,
        'f': lambda x, y: 30 * y * (x - 0.2) * (x - 0.7),
        'df_x': lambda x, y: y * (60 * x - 27),
        'df_xx': lambda x, y: y * 60,
        'df_y': lambda x, y: 30 * (x - 0.7) * (x - 0.2),
        'df_yy': lambda x, y: 0,
        'df_xy': lambda x, y: 60 * x - 27,
        'original_func': lambda x: 0.1 * math.exp(x * (10 * x * x - 13.5 * x + 4.2)),
    },
    '3': {
        'y0': 0.5,
        'f': lambda x, y: -20 * y * y * (x - 0.4),
        'df_x': lambda x, y: -20 * y * y,
        'df_xx': lambda x, y: 0,
        'df_y': lambda x, y: -20 * 2 * y * (x - 0.4),
        'df_yy': lambda x, y: -20 * 2 * (x - 0.4),
        'df_xy': lambda x, y: -20 * 2 * y,
        'original_func': lambda x: 0.1 / (x * x - 0.8 * x + 0.2),
    }
}
task = tasks[str(GROUP)]

methods = [
    (explicit_euler, 'r8-'),
    (implicit_euler, 'r--'),
    (euler_with_recount, 'rs-'),
    (cauchy, 'b8-'),
    (runge_kutta, 'yd-'),
    (extrapolation_adams, 'ko-'),
    (taylor_3, 'g^-'),
    (taylor_4, 'gs-'),
]


def draw_all(n):
    plt.close()
    x_points = np.linspace(0, 1, n)
    for method, color in methods:
        plt.plot(x_points, method(x_points, 1 / (n+1), **task), color, label=method.__doc__)
    plt.plot(x_points, list(map(task['original_func'], x_points)), 'b*-', label='Original')

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

l1 = Label(frame, text="Лабораторная работа №1.\n3 вариант.",
           justify=CENTER, font=("Helvetica", 12), bd=10, bg=color1)
l1.pack()

l = Label(frame, text="""f(x, y) = -20 * y^2 * (x - 0.4),\ty0 = 0.5
Решение численными методами
    1) Эйлера (явный)
    2) Эйлера (неявный)
    3) Эйлера (с пересчётом)
    4) Коши
    5) Рунге-Кутта 4 порядка
    6) Экстраполяционный метода Адамса (k=2)
    7) Тейлора (3 порядка)
    8) Тейлора (4 порядка),
а также точное решение.""", justify=LEFT, font=("Helvetica", 12), bd=10, bg=color2)
l.pack()

l2 = Label(frame, text="\nВыберите количество точек:", justify=CENTER, font=("Helvetica", 12), bd=0, bg=color1)
l2.pack()

w = Scale(frame, from_=20, to=500, resolution=10, length=300, bg=color1, borderwidth=0,
          relief=GROOVE, orient=HORIZONTAL, highlightthickness=0)
w.pack()

Label(frame, text='\n', bg=color1).pack()

button = Button(frame, text="Нарисовать график", font=("Helvetica", 12), bg=color2, command=lambda: draw_all(int(w.get())))
button.pack()

window.mainloop()
