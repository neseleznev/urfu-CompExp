#! /usr/bin/python3
# Nikita Seleznev
import math

def find_segment(f, eps):
    a = 0
    b = eps
    
    while (f(a) * f(b) > 0):
        a *= 2
        b *= 2
        
    return a, b

def dichotomy(f, a, b, eps):
    """Dichotomy:     """

    x = []
    while True:
        mid = (a + b) / 2.0
        x.append(mid)

        if abs(a - b) < eps:
            return x

        if f(a) * f(mid) < 0:
            b = mid
        elif f(mid) * f(b) < 0:
            a = mid


def fixed_chords(f, f_second_derivative, a, b, eps, m):
    """Fixed chords:   """

    x = []
    if f_second_derivative((a + b) / 2) > 0:
        x.append(a)  # x[0]
        x.append(b)  # x[i]
    else:
        x.append(b)  # x[0]
        x.append(a)  # x[i]
    i = 1

    while abs(f(x[i])/m) > eps:
        # x[i+1] =
        x.append( x[i] - (f(x[i]) * (x[i] - x[0])) / (f(x[i]) - f(x[0])) )
        i += 1
    return x[2:]


def floating_chords(f, f_second_derivative, a, b, eps, m):
    """Floating chords: """

    x = []
    if f_second_derivative((a + b) / 2) > 0:
        x.append(a)  # x[0]
        x.append(b)  # x[i]
    else:
        x.append(b)  # x[0]
        x.append(a)  # x[i]
    i = 1

    while abs(f(x[i])/m) > eps:
        # x[i+1] =
        x.append( x[i] - (f(x[i]) * (x[i] - x[i-1])) / (f(x[i]) - f(x[i-1])) )
        i += 1
    return x[2:]


def parabols(f, f_derivative, a, b, eps, M, m):
    x = []
    x.append(a)
    x.append(b)
       
        
    while ( abs(((x[-1] - x[-2]) **2) * M/m) > eps):
    
        c = (a + b) / 2;
        fa = f(a)
        fb = f(b)
        fc = f(c);

        A = ((fb - fc) / (b - c) - (fc - fa) / (c - a)) / (b - a);
        B = (fc - fa) / (c - a) + A* (a - c);
        C = f(a);

        x1 = a - (2 * C) / (B + (B * B - 4 * A * C) ** 0.5);
        x2 = a - (2 * C) / (B - (B * B - 4 * A * C) ** 0.5);
       

        if (x1 >= a and x1 <= b):
            x.append(x1)
        else:
            x.append(x2)

        if ( f(a) * f(x[-1]) < 0 ):
                b = x[-1]
        else:
                a = x[-1]
           
    
    return x[2:]
    
def newton(f, f_derivative, a, b, eps, M, m):
    """Newton method: """

    x = []
    x.append(a)  # x[0]
    x.append(b)  # x[i]
    i = 1

    while True:
        # x[i+1] =
        x.append( x[i] - f(x[i])/f_derivative(x[i]) )
        i += 1
        if abs(((x[i] - x[i-1]) ** 2) * M/(2*m)) < eps:
            break
    return x[2:]


def trim(floats, significant_position):
    p = significant_position
    return [(value // (10**-p)) / (10**p) for value in floats]


def print_results(method, f, a, b, eps, digits):
    values, root, iterations = method(f, a, b, eps)
    print("%s \t %.*f \t %d iterations." % (
        method.__doc__, digits, root, iterations))
    print("Values: " + str(trim(values, digits + 2)) + "\n\n")


def print_results_2(method, f, f_derivative, a, b, eps, digits):
    values, root, iterations = method(f, f_derivative, a, b, eps)
    print("%s \t %.*f \t %d iterations." % (
        method.__doc__, digits, root, iterations))
    print("Values: " + str(trim(values, digits + 2)) + "\n")


def main():
    GROUP = 4
    DIGITS = 8  # precision
    EPS = 0.5 * 1e-5
    print(EPS)
    group = [
        [
            lambda x: 1 + math.sin(x) - 1.2 * math.exp(-x),  # f(x)
            lambda x: math.cos(x) + 1.2 * math.exp(-x),      # f'(x)
            lambda x: -math.sin(x) - 1.2 * math.exp(-x),     # f''(x)
        ],
        [
            lambda x: math.log(x, 10) - 0.19/x,
            lambda x: (x+0.437491)/((x ** 2) * math.log(10, 2.718281828459) ),
            lambda x: (-0.434294*x - 0.38)/(x ** 3),
        ],
        [
            lambda x: x,
            lambda x: x,
            lambda x: x,
        ],
        [
            lambda x: 2*math.cos(x) - math.exp(x),
            lambda x: -2*math.sin(x) - math.exp(x),
            lambda x: -2*math.cos(x) - math.exp(x),
        ],
        [
            lambda x: x,
            lambda x: x,
            lambda x: x,
        ]
    ]

    print("\n\n\n")

    func = group[GROUP-1][0]
    first_dirivative = group[GROUP-1][1]
    second_dirivative = group[GROUP-1][2]
    
    a, b = find_segment(group[GROUP-1][0], EPS)
    
    print("start segment:")
    print(a, ":", b)
    print()
    
    m1 = min(group[GROUP-1][1](a), group[GROUP-1][1](b))
    M1 = max(group[GROUP-1][1](a), group[GROUP-1][1](b))
    M2 = max(group[GROUP-1][2](a), group[GROUP-1][2](b))
    
    #m1 = 0.373974
    #M2 = 0.814294
    
    print("M2 & m1:")
    print(M2, m1)
    
    real_root = 0.09327104616262699
    
    results = dichotomy(func, a, b, EPS)
    # results = [res - real_root for res in results]
    print("dichotomy:")
    print("num of iterations: ", len(results))
    print(trim(( results ), DIGITS))
    print()
    
    results = fixed_chords(func, second_dirivative, a, b, EPS, m1)
    # results = [res - real_root for res in results]
    print("fixed_chords:")
    print("num of iterations: ", len(results))
    print(trim(( results ), DIGITS))
    print()
    
    results = floating_chords(func, second_dirivative, a, b, EPS, m1)
    # results = [res - real_root for res in results]
    print("floating_chords:")
    print("num of iterations: ", len(results))
    print(trim(( results ), DIGITS))
    print()
    
    results = newton(func, first_dirivative, a, b, EPS, M2, m1)
    # results = [res - real_root for res in results]
    print("newton:")
    print("num of iterations: ", len(results))
    print(trim(( results ), DIGITS))
    print()
    
    results = parabols(func, first_dirivative, a, b, EPS, M2, m1)
    # results = [res - real_root for res in results]
    print("parabols:")
    print("num of iterations: ", len(results))
    print(trim(( results ), DIGITS))
    print()
    
    #print_results(dichotomy, group[GROUP-1][0], a, b, EPS, DIGITS)
    #
    #for method in [fixed_chords, floating_chords]:
    #    print_results_2(method, group[GROUP-1][0], group[GROUP-1][2], 0, 1, EPS, DIGITS)
    #
    #print_results_2(newton, group[GROUP-1][0], group[GROUP-1][1], 0, 1, EPS, DIGITS)
    #print_results_2(parabol, group[GROUP-1][0], group[GROUP-1][1], 0, 1, EPS, DIGITS)


if __name__ == '__main__':
    main()