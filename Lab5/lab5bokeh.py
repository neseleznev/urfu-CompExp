# Anna Mokrushina & RomanDubinin drawer
from bokeh.plotting import figure, output_file, show
from math import e
import numpy as np

from lab5 import L, f, real_func
from methods import runge_kutta, explicit_euler, recount_euler, get_next_m, progon


def main():

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
        y_points = [methods[i]["f"](x_points, h, y_0, last_m, f)]
        next_m = get_next_m(x_points, f, methods[i]["f"], y_0, right_boundary, last_m, h)

        while abs(next_m - last_m) > 0.01:
            y_points.append(methods[i]["f"](x_points, h, y_0, next_m, f))
            last_m = next_m
            next_m = get_next_m(x_points, f, methods[i]["f"], y_0, right_boundary, last_m, h)

        p.line(x=x_points, y=y_points[-1], color=colors[i], legend=methods[i]["name"])

    output_file("{0}_points.html".format(n))
    p.legend.orientation = "top_left"
    show(p)

if __name__ == '__main__':
    main()
