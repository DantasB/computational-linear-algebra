import math

from src.utils.matrix_utils import step_forward_derivate, step_backward_derivate, central_derivate


def forward_richard_extrapolation_derivate(function, x, delta, p):
    d1 = step_forward_derivate(function, x, delta)
    second_delta = delta/2
    d2 = step_forward_derivate(function, x, second_delta)
    q = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Forward:  " + str(result))


def backward_richard_extrapolation_derivate(function, x, delta, p):
    d1 = step_backward_derivate(function, x, delta)
    second_delta = delta/2
    d2 = step_backward_derivate(function, x, second_delta)
    q = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Backward: " + str(result))


def central_richard_extrapolation_derivate(function, x, delta, p):
    d1 = central_derivate(function, x, delta)
    second_delta = delta/2
    d2 = central_derivate(function, x, second_delta)
    q = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Central:  " + str(result))
