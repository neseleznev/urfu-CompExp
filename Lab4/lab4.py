from bokeh.plotting import figure, output_file, show
import bokeh
import itertools
import numpy as np
import math
import sys

Y0 = 0.1


def function_to_integrate(x, y):
    return 30 * y * (x - 0.2) * (x - 0.7)


def df_x(x, y):
    return y * (60 * x - 27)


def df_xx(x, y):
    return y * 60


def df_y(x, y):
    return 30 * (x - 0.7) * (x - 0.2)


def df_yy(x, y):
    return 0


def df_xy(x, y):
    return 60 * x - 27


def original_func(x):
    return 0.1 * math.exp(x * (10 * x * x - 13.5 * x + 4.2))


def explicit_euler(f, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        next_y = ys[k] + f(xs[k], ys[k]) * h
        ys.append(next_y)

    return ys[:-1]


def implicit_euler(f, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h
        next_y = ys[k] + f(xs[k], subsidiary_y) * h
        ys.append(next_y)

    return ys[:-1]


def cauchy(f, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h / 2
        next_y = ys[k] + f(xs[k] + h / 2, subsidiary_y) * h
        ys.append(next_y)

    return ys[:-1]


def euler_with_recount(f, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h
        next_y = ys[k] + (f(xs[k], ys[k]) + f(xs[k], subsidiary_y)) * h / 2
        ys.append(next_y)

    return ys[:-1]


def runge_kutta(f, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        k1 = h * f(xs[k], ys[k])
        k2 = h * f(xs[k] + h / 2, ys[k] + k1 / 2)
        k3 = h * f(xs[k] + h / 2, ys[k] + k2 / 2)
        k4 = h * f(xs[k] + h, ys[k] + k3)

        next_y = ys[k] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        ys.append(next_y)

    return ys[:-1]


def extrapolation_adams(f, xs, y0, h):
    ys = [y0, y0 + h * f(xs[0], y0)]
    for k in range(1, len(xs)):
        next_y = ys[k] + (
             3.0/2.0 * f(xs[k], ys[k]) - 1.0/2.0 * f(xs[k - 1], ys[k - 1])
        ) * h
        ys.append(next_y)

    return ys[:-1]


def taylor_3(f, df_x, df_y, xs, y0, h):
    ys = [y0]
    for k in range(0, len(xs)):
        next_y = ys[k] + f(xs[k], ys[k])*h + (
             df_x(xs[k], ys[k]) + df_y(xs[k], ys[k])*f(xs[k], ys[k])
        ) * (h ** 2) / 2
        ys.append(next_y)

    return ys[:-1]


def taylor_4(f, df_x, df_y, df_xx, df_yy, df_xy, xs, y0, h):
    ys = [y0]
    for k in range(len(xs)):
        second_summand = f(xs[k], ys[k]) * h
        third_summand = (
            df_x(xs[k], ys[k]) + df_y(xs[k], ys[k]) * f(xs[k], ys[k])
        ) * (h ** 2) / 2
        fourth_summand = (
             df_xx(xs[k], ys[k]) + 2 * f(xs[k], ys[k]) * df_xy(xs[k], ys[k]) +
             df_yy(xs[k], ys[k]) * (f(xs[k], ys[k]) ** 2) +
             df_y(xs[k], ys[k]) * (
                 df_x(xs[k], ys[k]) + df_y(xs[k], ys[k]) * f(xs[k], ys[k]))
        ) * (h ** 3) / 3
        next_y = ys[k] + second_summand + third_summand + fourth_summand
        ys.append(next_y)

    return ys[:-1]


n = int(sys.argv[1])

output_file("{0}_points.html".format(n))

p = figure(title="30*y*(x - 0.2)*(x - 0.7), " + str(n) + "points", plot_width=1100, plot_height=600)

x_points = np.linspace(0, 1, n, endpoint=True)

y_points = explicit_euler(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="red", legend="Эйлер")
p.circle(x=x_points, y=y_points, color="red", legend="Эйлер")

y_points = cauchy(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="blue", legend="Коши")

y_points = implicit_euler(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="red", line_dash=[4, 4], legend="неявный Эйлер")

y_points = euler_with_recount(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="red", legend="Эйлер с пересчётом")
p.square(x=x_points, y=y_points, color="red", legend="Эйлер с пересчётом")

y_points = runge_kutta(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="purple", legend="Рунге-Кутты 4 го порядка")

y_points = extrapolation_adams(function_to_integrate, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="orange", legend="экстраполяционный метод Адамса(k=2)")

y_points = taylor_3(function_to_integrate, df_x, df_y, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="green", legend="Тейлор 3 го порядка")
p.triangle(x=x_points, y=y_points, color="green", legend="Тейлор 3 го порядка")

y_points = taylor_4(function_to_integrate, df_x, df_y, df_xx, df_yy, df_xy, x_points, Y0, 1 / n)
p.line(x=x_points, y=y_points, color="green", legend="Тейлор 4 го порядка")
p.square(x=x_points, y=y_points, color="green", legend="Тейлор 4 го порядка")

y_points = [original_func(x) for x in x_points]
p.line(x=x_points, y=y_points, color="black", legend="origin")
p.circle(x=x_points, y=y_points, color="black", legend="origin")

p.legend.orientation = "top_left"
show(p)
