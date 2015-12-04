#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2015

__author__ = 'Nikita'


class MethodError(Exception):
    pass


def trim(value, accuracy):
    return int((value * (10 ** accuracy))) / (10 ** accuracy)


def get_indexes_of_no_zero_items(list_, accuracy):
    return set(i for i, element in enumerate(list_)
               if abs(trim(element, accuracy - 1)) > 0)


def matrix_to_str(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    return '\n'.join(table)


def multiply_matrices(A, B, accuracy):
    if len(A[0]) != len(B):
        raise ValueError("Cannot multiply matrices. Incorrect dimensions.")

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(len(B[0]))] for col in range(len(A))]

    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(A[0])):
                C[i][j] += trim(A[i][k] * B[k][j], accuracy)
    return C


def swap_matrix_lines(matrix, string_num_1, string_num_2):
    temp = matrix[string_num_1]
    matrix[string_num_1] = matrix[string_num_2]
    matrix[string_num_2] = temp


def find_max_cell(matrix, start_string_num):
    max_val = matrix[start_string_num][0]
    res_i = start_string_num
    res_j = 0
    for i in range(start_string_num, len(matrix)):
        for j in range(0, len(matrix[i])):
            if abs(matrix[i][j]) > max_val:
                max_val = matrix[i][j]
                res_i = i
                res_j = j
    return res_i, res_j


def distance_R3_vectors(m1, m2):
    return (
        (m1[0][0] - m2[0][0]) ** 2 +
        (m1[1][0] - m2[1][0]) ** 2 +
        (m1[2][0] - m2[2][0]) ** 2
    ) ** (1/2)
