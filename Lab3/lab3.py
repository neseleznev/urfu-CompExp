#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nikita Seleznev, 2015
from gauss import compact_gauss_scheme, main_item_gauss_scheme
from gauss.utils import MethodError, matrix_to_str, trim, distance_R3_vectors

N = 14
M = 20 - N
A = 0.1 * M + 0.01 * N
B = 0.2 * M + 0.02 * N + 0.01*M*N + 0.001*N*N

set_of_equations = [
    [1.2345, 3.1415, 1,         7.5175 + A],
    [2.3456, 5.9690, 0,         14.2836],
    [3.4567, 2.1828, 2 + 0.1*N, 7.8223 + B]
]

X_asterisk = [
    [1],
    [2],
    [0.1*M + 0.01*N],
]


def do_method_investigation(method, accuracy_list, verbose=True):
    for accuracy in accuracy_list:
        copy_of_system = []

        for i in range(len(set_of_equations)):
            copy_of_system.append([])
            for j in range(len(set_of_equations[i])):
                copy_of_system[i].append(
                    trim(set_of_equations[i][j], accuracy))

        try:
            print("\nAccuracy: ", accuracy)
            x_by_method = method(copy_of_system, accuracy, verbose)
            print("Answer:\n" + matrix_to_str(x_by_method) + "\n")
            print("Distance R3: " + str(distance_R3_vectors(x_by_method, X_asterisk)))
        except MethodError as e:
            print(e)


def main():
    import subprocess
    subprocess.call("clear", shell=True)
    subprocess.call("cls", shell=True)

    print("\n[I]\t\t\tReal answer(X*):\n")
    print(matrix_to_str(X_asterisk))
    print("Proof:")
    import sys
    sys.stdout.write('http://www.wolframalpha.com/input/?i=gauss [')
    for row in set_of_equations[:-1]:
        sys.stdout.write(str(row) + ", ")
    sys.stdout.write(str(set_of_equations[-1]))
    sys.stdout.write(']\n\n')

    print("=" * 80)

    print("[II]\t\t\tCompact Gauss Scheme: ")
    do_method_investigation(compact_gauss_scheme, accuracy_list=[2, 4, 6])
    print("=" * 80)

    print("[III]\t\t\tGauss Scheme with main element: ")
    do_method_investigation(main_item_gauss_scheme, accuracy_list=[2, 4, 6])


if __name__ == '__main__':
    main()
