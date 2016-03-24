#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
import sys
from bokeh.plotting import figure, output_file, show
import inspect
import numpy as np
import math

from Lab4.methods import explicit_euler, euler_with_recount, cauchy, \
    runge_kutta, implicit_euler, taylor_3, taylor_4, extrapolation_adams

GROUP = 4

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
    '4': {
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

n = int(sys.argv[1])
output_file("{0}_points.html".format(n))


x_points = np.linspace(0, 1, n, endpoint=True)
Y0, f, df_x, df_xx, df_y, df_yy, df_xy, solution = task['y0'],\
    task['f'], task['df_x'], task['df_xx'], task['df_y'],\
    task['df_yy'], task['df_xy'], task['original_func']

p = figure(
    title='{} {} points'.format(
        inspect.getsource(f).split(':')[-1].strip(), n),
    plot_width=1100,
    plot_height=600
)

[draw(x=x_points, y=explicit_euler(f, x_points, Y0, 1 / n),
      color="red", legend="Эйлер"
      ) for draw in [p.line, p.circle]]

p.line(x=x_points, y=cauchy(f, x_points, Y0, 1 / n),
       color="blue", legend="Коши")

p.line(x=x_points, y=implicit_euler(f, x_points, Y0, 1 / n),
       color="red", legend="Неявный метод Эйлера")

[draw(x=x_points, y=euler_with_recount(f, x_points, Y0, 1 / n),
      color="red", legend="Эйлера с пересчётом"
      ) for draw in [p.line, p.square]]

p.line(x=x_points, y=runge_kutta(f, x_points, Y0, 1 / n),
       color="purple", legend="Рунге-Кутты 4 го порядка")

p.line(x=x_points, y=extrapolation_adams(f, x_points, Y0, 1 / n),
       color="orange", legend="Экстраполяционный метод Адамса(k=2)")

[draw(x=x_points, y=taylor_3(f, df_x, df_y, x_points, Y0, 1 / n),
      color="green", legend="Тейлор 3 го порядка"
      ) for draw in [p.line, p.triangle]]

[draw(x=x_points, y=taylor_4(f, df_x, df_y, df_xx, df_yy, df_xy, x_points, Y0, 1 / n),
      color="green", legend="Тейлор 4 го порядка"
      ) for draw in [p.line, p.square]]

[draw(x=x_points, y=list(map(solution, x_points)),
      color="black", legend="origin"
      ) for draw in [p.line, p.circle]]

p.legend.location = "top_left"
show(p)
