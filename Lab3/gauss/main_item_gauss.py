#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2015
import copy
from .utils import matrix_to_str, swap_matrix_lines, find_max_cell, \
    MethodError, get_indexes_of_no_zero_items, trim


def solve_any_transformed_system(set_of_equations, accuracy):
    x_column = []
    for _ in range(len(set_of_equations[0]) - 1):
        x_column.append([0])

    computed_elements_of_x = []
    num_of_computed_x_on_last_iteartion = 0
    no_zero_indexes_on_last_iteration = []

    while len(computed_elements_of_x) != len(x_column):
        for equation in set_of_equations:
            no_zero_elements_in_current_row = \
                get_indexes_of_no_zero_items(equation[:-1], accuracy)

            new_no_zero_elements = list(
                no_zero_elements_in_current_row -
                set(no_zero_indexes_on_last_iteration)
            )
            if len(new_no_zero_elements) != 1:
                continue

            j = new_no_zero_elements[0]

            sum_ = 0
            for k in range(0, j):
                sum_ += trim(equation[k] * x_column[k][0], accuracy)
            for k in range(j + 1, len(set_of_equations[0]) - 1):
                sum_ += trim(equation[k] * x_column[k][0], accuracy)

            x_column[j][0] = trim((equation[-1] - sum_) / equation[j],
                                  accuracy)

            computed_elements_of_x.append(j)
            num_of_computed_x_on_last_iteartion += 1
            no_zero_indexes_on_last_iteration = no_zero_elements_in_current_row
            break

        if (num_of_computed_x_on_last_iteartion ==
           len(computed_elements_of_x) - 1):
            raise MethodError("Method is unable to find solution!")
    return x_column


def do_forward_step_of_gauss_by_main_element(matrix, accuracy):
    set_of_equations = copy.deepcopy(matrix)

    for start_string_index in range(len(set_of_equations)):
        A = [line[:-1] for line in set_of_equations]

        main_i, main_j = find_max_cell(A, start_string_index)
        swap_matrix_lines(set_of_equations, start_string_index, main_i)

        for equation in set_of_equations[start_string_index + 1:]:
            multiplicator = trim(-equation[main_j] /
                                 set_of_equations[start_string_index][main_j],
                                 accuracy)
            for j in range(len(equation)):
                equation[j] = trim(
                    equation[j] +
                    set_of_equations[start_string_index][j]*multiplicator,
                    accuracy
                )
    return set_of_equations


def main_item_gauss_scheme(matrix, accuracy, verbose=True):
    set_of_equations = do_forward_step_of_gauss_by_main_element(
        matrix, accuracy)
    if verbose:
        print("Set of equations after forward step of method:\n" +
              matrix_to_str(set_of_equations) + "\n")
    return solve_any_transformed_system(set_of_equations, accuracy)
