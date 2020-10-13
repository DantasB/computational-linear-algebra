from importlib.machinery import SourceFileLoader

import math

Matrix_Utils = SourceFileLoader("Matrix_Utils", "/home/bdantas/Área de Trabalho/ALC_Lists/Utils/Matrix_Utils.py").load_module()


def forward_richard_extrapolation_derivate(function, x, delta, p):
    d1           = Matrix_Utils.step_forward_derivate(function, x, delta)
    second_delta = delta/2
    d2           = Matrix_Utils.step_forward_derivate(function, x, second_delta)
    q            = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Forward:  " + str(result))


def backward_richard_extrapolation_derivate(function, x, delta, p):
    d1           = Matrix_Utils.step_backward_derivate(function, x, delta)
    second_delta = delta/2
    d2           = Matrix_Utils.step_backward_derivate(function, x, second_delta)
    q            = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Backward: " + str(result))


def central_richard_extrapolation_derivate(function, x, delta, p):
    d1           = Matrix_Utils.central_derivate(function, x, delta)
    second_delta = delta/2
    d2           = Matrix_Utils.central_derivate(function, x, second_delta)
    q            = delta/second_delta

    result = d1 + (d1 - d2)/(q**(-p)-1)

    print("Richard Extrapolation - Central:  " + str(result))
