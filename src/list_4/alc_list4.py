import math
from src.utils.matrix_utils import get_transposed_matrix, broyden_method_helper, sum_matrixes, vector_multiplication, multiply_matrixes, multiply_matrix_scalar, derivate, inverse_interpolation_helper, get_jacobian_matrix, get_f_vector, get_inverse_matrix, multiply_matrix_vector, sum_vectors, vector_norm, get_inverse_matrix, get_f_vector

STEPS = 100
TOL = 10**-4


def bissection_method(function, a, b):
    if function(a)*function(b) > 0:
        print("Error")
        return

    c = a
    while((b - a) >= TOL):

        c = (a+b)/2
        if(function(c) == 0):
            break

        if (function(c)*function(a) < 0):
            b = c
        else:
            a = c

    print("Root: " + str(c))


def newton_method(function, x0):
    x = x0

    for _ in range(STEPS):
        xi = x - function(x) / derivate(function, x)
        residue = math.fabs(xi - x)

        if(residue < TOL):
            print("Root: " + str(xi))
            return

        x = xi

    print("Convergence not reached")


def secant_method(function, x0):
    delta_x = 0.001
    first_x = x0
    actual_x = first_x + delta_x
    first_function = function(first_x)

    for _ in range(STEPS):
        actual_function = function(actual_x)
        next_x = actual_x - \
            float((actual_function*(actual_x - first_x)) /
                  (actual_function - first_function))
        residue = abs(next_x - actual_x)

        if (residue < TOL):
            print("Root: " + str(actual_x))
            return

        first_function = actual_function
        first_x = actual_x
        actual_x = next_x

    print("Convergence not reached")


def inverse_interpolation_method(function, x1, x2, x3):
    x0 = 10**(36)
    x = [x1, x2, x3]
    y = [function(x1), function(x2), function(x3)]

    for _ in range(STEPS):
        xi = inverse_interpolation_helper(
            x[0], x[1], x[2], y[0], y[1], y[2])
        residue = abs(xi - x0)

        if (residue < TOL):
            print("Root: " + str(xi))
            return

        x0 = xi

        y_max = max(y)
        if(y_max == y[0]):
            x1 = x0
            x = [x1, x2, x3]
            x.sort()
            y = [function(x[0]), function(x[1]), function(x[2])]
            continue

        elif(y_max == y[0]):
            x2 = x0
            x = [x1, x2, x3]
            x.sort()
            y = [function(x[0]), function(x[1]), function(x[2])]
            continue

        else:
            x3 = x0
            x = [x1, x2, x3]
            x.sort()
            y = [function(x[0]), function(x[1]), function(x[2])]
            continue

    print("Convergence not reached")


def non_linear_newton_method(functions_list, first_solution):
    x = first_solution

    for _ in range(STEPS):
        jacobian_matrix = get_jacobian_matrix(functions_list, x)
        vector_f = get_f_vector(functions_list, x)
        jacobian_inverse = get_inverse_matrix(jacobian_matrix)
        if(jacobian_inverse == 0):
            print("Error")
            break

        delta_x = multiply_matrix_vector(
            jacobian_inverse, vector_f)
        delta_x = [x*(-1) for x in delta_x]
        x = sum_vectors(x, delta_x)

        residue = vector_norm(delta_x)/vector_norm(x)

        if (residue < TOL):
            print("Solution: " + str(x))
            return

    print("Convergence not reached")


def broyden_method(functions_list, first_solution):
    dimension = len(first_solution)
    x = first_solution
    jacobian_b = get_jacobian_matrix(functions_list, x)
    jacobian_a = [[float(i == j) for j in range(dimension)]
                  for i in range(dimension)]

    for _ in range(STEPS):

        for k in range(dimension):
            for j in range(dimension):
                jacobian_a[k][j] = jacobian_b[k][j]

        vector_f = get_f_vector(functions_list, x)
        jacobian_inverse = get_inverse_matrix(jacobian_a)
        if(jacobian_inverse == 0):
            print("Error")
            break

        delta_x = multiply_matrix_vector(
            jacobian_inverse, vector_f)
        delta_x = [x*(-1) for x in delta_x]

        x = sum_vectors(x, delta_x)

        vector_f2 = get_f_vector(functions_list, x)

        yK = sum_vectors(vector_f2, [x*(-1) for x in vector_f])
        residue = vector_norm(delta_x)/vector_norm(x)
        if (residue < TOL):
            print("Solution: " + str(x))
            return

        product = multiply_matrix_vector(jacobian_b, delta_x)
        product = [x*(-1) for x in product]

        denominator = vector_multiplication(delta_x, delta_x)
        top = sum_vectors(yK, product)
        numerator = broyden_method_helper(top, delta_x)

        for k in range(len(jacobian_b)):
            for j in range(len(jacobian_b)):
                jacobian_b[k][j] = numerator[k][j]/denominator

        jacobian_b = sum_matrixes(jacobian_a, jacobian_b)

    print("Convergence not reached")


def non_linear_mmq(functions_list, first_solution):
    x = first_solution

    for _ in range(STEPS):
        jacobian = get_jacobian_matrix(functions_list, x)
        transposed_jacobian = get_transposed_matrix(jacobian)
        vector_f = get_f_vector(functions_list, x)

        a = get_inverse_matrix(
            multiply_matrixes(transposed_jacobian, jacobian))
        if(a == 0):
            print("Error")
            break

        b = multiply_matrix_scalar(a, -1)
        c = multiply_matrix_vector(transposed_jacobian, vector_f)
        delta_b = multiply_matrix_vector(b, c)
        x = sum_vectors(x, delta_b)

        residue = vector_norm(delta_b)/vector_norm(x)

        if (residue < TOL):
            print("Solution: " + str(x))
            return

    print("Convergence not reached")
