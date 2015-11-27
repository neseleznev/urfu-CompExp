#!/usr/bin/env python3
# -*- coding: utf8 -*-
# Nikita Seleznev, 2015
import math


def dichotomy(functions, m_and_Ms, a, b, eps):
    """Dichotomy:     """

    f = functions[0]
    x = []
    while True:
        mid = (a + b) / 2.0
        x.append(mid)

        if abs(a - b) < eps:
            return x, mid, len(x)

        if f(a) * f(mid) < 1e-15:
            b = mid
        elif f(mid) * f(b) < 1e-15:
            a = mid


def fixed_chords(functions, m_and_Ms, a, b, eps):
    """Fixed chords:   """

    f, f_sec_derivative = functions
    m = m_and_Ms[0]
    x = []
    if f(a) * f_sec_derivative(a) > 1e-15:
        x.append(a)  # x[0]
        x.append(b)  # x[1]
    elif f(b) * f_sec_derivative(b) > 1e-15:
        x.append(b)  # x[0]
        x.append(a)  # x[1]
    else:
        raise ValueError("f must retain it's second derivative sign for [a;b]")
    i = 1

    while abs(f(x[i]) / m) > eps:  # abs(x[i] - x[i-1]) > eps:
        # x[i+1] =
        x.append(x[i] - (f(x[i])*(x[i] - x[0])) / (f(x[i])-f(x[0])))
        i += 1
    return x[2:], x[i], i-1


def floating_chords(functions, m_and_Ms, a, b, eps):
    """Floating chords: """

    f, f_sec_derivative = functions
    m = m_and_Ms[0]
    x = []
    if f(a) * f_sec_derivative(a) > 1e-15:
        x.append(a)  # x[0]
        x.append(b)  # x[1]
    elif f(b) * f_sec_derivative(b) > 1e-15:
        x.append(b)  # x[0]
        x.append(a)  # x[1]
    else:
        raise ValueError("f must retain it's second derivative sign for [a;b]")
    i = 1

    while abs(f(x[i]) / m) > eps:  # abs(x[i] - x[i-1]) > eps:
        # x[i+1] =
        x.append(x[i] - (f(x[i]) * (x[i] - x[i-1])) / (f(x[i]) - f(x[i-1])))
        i += 1
    return x[2:], x[i], i-1


def muller(functions, m_and_Ms, a, b, eps):
    """Muller method:   """

    f = functions[0]
    m, M = m_and_Ms
    c = 0.0
    x0, x = a, b
    i = 1

    fa, fb, fc = 0.0, 0.0, 0.0

    ans = [x]

    while abs(((x - x0) ** 2) * M/m) > eps:  # abs(x - x0) > eps:
        x0 = x
        c = (a + b) / 2.0
        fa = f(a)
        fb = f(b)
        fc = f(c)

        A = ((fb - fc) / (b - c) - (fc - fa) / (c - a)) / (b - a)
        B = (fc - fa) / (c - a) + A*(a - c)
        C = f(a)

        x1 = a - (2*C) / (B + math.sqrt(B*B - 4*A*C))
        x2 = a - (2*C) / (B - math.sqrt(B*B - 4*A*C))

        if a <= x1 <= b:
            x = x1
        else:
            x = x2
        ans.append(x)

        if f(a)*f(x) < 0:
            b = x
        else:
            a = x
        i += 1
    return ans[1:], ans[-1], i-1


def newton(functions, m_and_Ms, a, b, eps):
    """Newton method: """

    f, f_derivative, f_sec_derivative = functions
    m, M = m_and_Ms

    x = []
    if f(a) * f_sec_derivative(a) > 1e-15:
        x.append(a)  # x[0]
    elif f(b) * f_sec_derivative(b) > 1e-15:
        x.append(b)  # x[0]
    else:
        raise ValueError("f must retain it's second derivative sign for [a;b]")
    x.append(x[0] - f(x[0]) / f_derivative(x[0]))  # x[1]
    i = 1

    while abs(((x[i] - x[i-1]) ** 2) * M/(2*m)) > eps:  # abs(x[i] - x[i-1]) > eps:
        # x[i+1] =
        x.append(x[i] - f(x[i]) / f_derivative(x[i]))
        i += 1
    return x[1:], x[i], i


def find_interval(f, eps):
    a, b = 0.0, eps
    while f(a) * f(b) > 0:
        b *= 2
    return a, b


def drange(start, stop, step):
    while start < stop:
        yield start
        start += step


def trim(floats, significant_position):
    p = significant_position
    return [(value // (10**-p)) / (10**p) for value in floats]


def print_results(method, functions, m_and_Ms, a, b, eps, digits):
    values, root, iterations = method(functions, m_and_Ms, a, b, eps)
    print("%s \t f(%.*f) = %.*f \t %d iterations." % (
        method.__doc__, digits, root, digits, functions[0](root), iterations))
    print("Values: " + str(trim(values, digits + 2)) + "\n")
    # print("F(x): " + str(trim([functions[0](v) for v in values], digits + 2)) + "\n")


def main():
    GROUP = 4
    DIGITS = 6  # precision
    EPS = 0.5 * 1e-5

    group = [
        [
            lambda x: 1 + math.sin(x) - 1.2 * math.exp(-x),  # f(x)
            lambda x: math.cos(x) + 1.2 * math.exp(-x),      # f'(x)
            lambda x: -math.sin(x) - 1.2 * math.exp(-x),     # f''(x)
        ],
        [
            lambda x: NotImplemented,
            lambda x: NotImplemented,
            lambda x: NotImplemented,
        ],
        [
            lambda x: NotImplemented,
            lambda x: NotImplemented,
            lambda x: NotImplemented,
        ],
        [
            lambda x: 2*math.cos(x) - math.exp(x),
            lambda x: -2*math.sin(x) - math.exp(x),
            lambda x: -2*math.cos(x) - math.exp(x),
        ],
        [
            lambda x: NotImplemented,
            lambda x: NotImplemented,
            lambda x: NotImplemented,
        ]
    ]

    print("\n\n\n")
    a, b = find_interval(group[GROUP-1][0], EPS)
    print("Interval [{}; {}]".format(a, b))


    m1 = min([group[GROUP-1][1](x) for x in drange(a, b, EPS)])
    M1 = max([group[GROUP-1][1](x) for x in drange(a, b, EPS)])
    m2 = min([group[GROUP-1][2](x) for x in drange(a, b, EPS)])
    M2 = max([group[GROUP-1][2](x) for x in drange(a, b, EPS)])

    print_results(method=dichotomy,
                  functions=(group[GROUP-1][0], ),
                  m_and_Ms=None,
                  a=a, b=b,
                  eps=EPS, digits=DIGITS)
    print_results(method=fixed_chords,
                  functions=(group[GROUP-1][0], group[GROUP-1][2]),
                  m_and_Ms=(m1, ),
                  a=a, b=b,
                  eps=EPS, digits=DIGITS)
    print_results(method=floating_chords,
                  functions=(group[GROUP-1][0], group[GROUP-1][2]),
                  m_and_Ms=(m1, ),
                  a=a, b=b,
                  eps=EPS, digits=DIGITS)
    print_results(method=muller,
                  functions=(group[GROUP-1][0], ),
                  m_and_Ms=(m1, M2, ),
                  a=a, b=b,
                  eps=EPS, digits=DIGITS)
    print_results(method=newton,
                  functions=(group[GROUP-1][0], group[GROUP-1][1], group[GROUP-1][1]),
                  m_and_Ms=(m1, M2, ),
                  a=a, b=b,
                  eps=EPS, digits=DIGITS)


if __name__ == '__main__':
    main()
