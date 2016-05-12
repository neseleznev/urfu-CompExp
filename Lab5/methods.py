#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2016
# Great thanks to github.com/RomanDubinin


def explicit_euler(x_points, h, y_0, z_0, func):
    """Explicit Euler"""
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points)-1):
        new_y = y_points[i] + h * z_points[i]
        new_z = z_points[i] + h * func(x_points[i], y_points[i], z_points[i])
        y_points.append(new_y)
        z_points.append(new_z)

    return y_points


def recount_euler(x_points, h, y_0, z_0, func):
    """Euler with recount"""
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points) - 1):
        new_y = y_points[i] + (h/2) * (2 * z_points[i] + h * func(x_points[i], y_points[i], z_points[i]))
        new_z = z_points[i] + (h / 2) * (
            func(x_points[i], y_points[i], z_points[i]) +
            func(
                x_points[i + 1],
                y_points[i] + h * z_points[i],
                z_points[i] + h * func(x_points[i], y_points[i], z_points[i])
            )
        )
        y_points.append(new_y)
        z_points.append(new_z)

    return y_points


def runge_kutta(x_points, h, y_0, z_0, func):
    """Runge-Kutta"""
    y_points = [y_0]
    z_points = [z_0]

    for i in range(len(x_points) - 1):
        k1 = h * func(x_points[i], y_points[i], z_points[i])
        k2 = h * func(x_points[i] + h / 2, y_points[i] + k1 / 2, z_points[i])
        k3 = h * func(x_points[i] + h / 2, y_points[i] + k2 / 2, z_points[i])
        k4 = h * func(x_points[i] + h, y_points[i] + k3, z_points[i])

        new_y = y_points[i] + h*z_points[i]
        new_z = z_points[i] + (k1 + 2*k2 + 2*k3 + k4)/6

        y_points.append(new_y)
        z_points.append(new_z)

    return y_points


def get_next_m(x_points, func, method, y_0, dy_n, last_m, h):
    y_points = method(x_points, h, y_0, last_m, func)
    fi = ((y_points[-1] - y_points[-2]) / h) - dy_n

    y_points = method(x_points, h, 0, 1, lambda x, y, d_y: y)
    d_fi = (y_points[-1] - y_points[-2]) / h

    return last_m - fi/d_fi


def progon(x_points, dy_n, h, l):
    lambda_ = [0]
    mu = [0]
    last_number = len(x_points) - 1
    for i in range(len(x_points)):
        lambda_.append(-1/(lambda_[i] - 2 - h*h))
        mu.append((h*h*(l*x_points[i]*(1-x_points[i]) + 2*l + 2) - mu[i]) / (lambda_[i] - 2 - h*h))

    y_points = [(h*dy_n + mu[last_number]) / (1 - lambda_[last_number])]*len(x_points)
    for i in range(last_number, 0, -1):
        y_points[i-1] = (lambda_[i-1]*y_points[i] + mu[i-1])

    return y_points
