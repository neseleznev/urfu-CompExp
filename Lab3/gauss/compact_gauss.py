#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2015
from .utils import MethodError, matrix_to_str, trim


def calculate_cell_by_gauss(i, j, A, B, C, accuracy):
    sum_ = 0
    for k in range(0, i):
        sum_ += trim(B[i][k] * C[k][j], accuracy)
    for k in range(i + 1, len(B[1])):
        sum_ += trim(B[i][k] * C[k][j], accuracy)

    if j > i:  # this element is in C
        C[i][j] = trim((A[i][j] - sum_) / B[i][i], accuracy)
    else:  # element in B
        B[i][j] = trim((A[i][j] - sum_) / C[i][i], accuracy)


def solve_up_triangle_system(matrix, b_column, accuracy):
    x_column = [[0] for _ in range(len(b_column))]

    for i in range(len(b_column))[::-1]:
        sum_ = 0
        for k in range(i):
            sum_ += trim(matrix[i][k] * x_column[k][0], accuracy)

        for k in range(i+1, len(b_column)):
            sum_ += trim(matrix[i][k] * x_column[k][0], accuracy)

        x_column[i][0] = trim((b_column[i][0] - sum_) /
                              matrix[i][i], accuracy)
    return x_column


def calculate_matrix_B_nad_C_by_gauss(A, accuracy):
    B = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ]
    C = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]
    ]

    height = len(A)
    weedth = len(A[0])

    try:
        for k in range(0, weedth):
            for i in range(k, height):
                calculate_cell_by_gauss(i, k, A, B, C, accuracy)
            for j in range(k + 1, weedth):
                calculate_cell_by_gauss(k, j, A, B, C, accuracy)
    except ZeroDivisionError:
        raise MethodError("Method can not find solution")

    return B, C


def compact_gauss_scheme(set_of_equations, accuracy, verbose=True):
    B, C = calculate_matrix_B_nad_C_by_gauss(set_of_equations, accuracy)

    y_column = [[row[-1]] for row in C]

    C = [line[:-1] for line in C]

    x_by_method = solve_up_triangle_system(C, y_column, accuracy)

    if verbose:
        print("Matrix B:\n" + matrix_to_str(B) + "\n")
        print("Matrix C:\n" + matrix_to_str(C) + "\n")
        print("column y:\n" + matrix_to_str(y_column) + "\n")
    return x_by_method
