#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
# Great thanks to github.com/RomanDubinin


def explicit_euler(xs, h, y0, f, **derivatives):
    """Explicit Euler"""
    ys = [y0]
    for k in range(len(xs)):
        next_y = ys[k] + f(xs[k], ys[k]) * h
        ys.append(next_y)

    return ys[:-1]


def implicit_euler(xs, h, y0, f, **derivatives):
    """Implicit Euler"""
    ys = [y0]
    for k in range(len(xs) - 1):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h
        next_y = ys[k] + f(xs[k + 1], subsidiary_y) * h
        ys.append(next_y)

    return ys


def cauchy(xs, h, y0, f, **derivatives):
    """Cauchy"""
    ys = [y0]
    for k in range(len(xs) - 1):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h / 2
        next_y = ys[k] + f(xs[k] + h / 2, subsidiary_y) * h
        ys.append(next_y)

    return ys


def euler_with_recount(xs, h, y0, f, **derivatives):
    """Euler with recount"""
    ys = [y0]
    for k in range(len(xs) - 1):
        subsidiary_y = ys[k] + f(xs[k], ys[k]) * h
        next_y = ys[k] + (f(xs[k], ys[k]) + f(xs[k], subsidiary_y)) * h / 2
        ys.append(next_y)

    return ys


def runge_kutta(xs, h, y0, f, **derivatives):
    """Runge-Kutta 4th"""
    ys = [y0]
    for k in range(len(xs)):
        k1 = h * f(xs[k], ys[k])
        k2 = h * f(xs[k] + h / 2, ys[k] + k1 / 2)
        k3 = h * f(xs[k] + h / 2, ys[k] + k2 / 2)
        k4 = h * f(xs[k] + h, ys[k] + k3)

        next_y = ys[k] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        ys.append(next_y)

    return ys[:-1]


def extrapolation_adams(xs, h, y0, f, **derivatives):
    """Extra Adams (k=2)"""
    ys = [y0, y0 + h * f(xs[0], y0)]
    for k in range(1, len(xs)):
        next_y = ys[k] + (
             3.0/2.0 * f(xs[k], ys[k]) - 1.0/2.0 * f(xs[k - 1], ys[k - 1])
        ) * h
        ys.append(next_y)

    return ys[:-1]


def taylor_3(xs, h, y0, f, df_x, df_y, **derivatives):
    """Taylor 3rd"""
    ys = [y0]
    for k in range(0, len(xs)):
        next_y = ys[k] + f(xs[k], ys[k])*h + (
             df_x(xs[k], ys[k]) + df_y(xs[k], ys[k])*f(xs[k], ys[k])
        ) * (h ** 2) / 2
        ys.append(next_y)

    return ys[:-1]


def taylor_4(xs, h, y0, f, df_x, df_y, df_xx, df_yy, df_xy, **derivatives):
    """Taylor 4rd"""
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
