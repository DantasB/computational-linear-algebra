import math

from src.utils import matrix_utils
from src.list_1 import alc_list1
from src.list_2 import alc_list2
from src.list_3 import alc_list3
from src.list_4 import alc_list4
from src.list_5 import alc_list5
from src.list_6 import alc_list6
from src.list_7 import alc_list7


def main(function, arg1, arg2=[], arg3=[], arg4=[], arg5=[], arg6=[]):
    # Função, menor valor, maior valor
    #alc_list4.bissection_method(function, arg1, arg2)

    # função, número inicial
    #alc_list4.newton_method(function, arg1)

    # função, número inicial
    #alc_list4.secant_method(function, arg1)

    # função, 3 números crescentes próximos à raiz
    #alc_list4.inverse_interpolation_method(function, arg1, arg2, arg3)

    # lista de funções, chute inicial de solução
    #alc_list4.non_linear_newton_method(function, arg1)

    # lista de funções, chute inicial de solução
    #alc_list4.broyden_method(function, arg1)

    # lista de funções reglin, vetor x, vetor y, chute inicial de solução
    alc_list4.non_linear_mmq(function, arg1)

    #função, a, b
    #alc_list5.estimate_integral_value(function, arg1, arg2)

    # função, a, b, sensibilidade, boleano (True = polinomial, False = gauss)
    #alc_list5.integrate(function, arg1, arg2, arg3, False)

    # função diferencial, t0, tf, delta = 0.1, condição inicial, escolha (0 = euler, 1 = RK2, Default = RK4)
    #alc_list6.first_order_edo_solver(function, arg1, arg2, arg3, arg4, arg5)

    # função diferencial de segunda ordem, t0, tf, delta = 0.1, condição inicial, condição inicial da derivada, escolha(true=Taylor, false=RKN)
    #alc_list6.second_order_edo_solver(function, arg1, arg2, arg3, arg4, arg5, arg6)

    #função, x, deltaX
    #print("Derivate              - Forward:  " + str(matrix_utils.step_forward_derivate(function, arg1, arg2)))

    #função, x, deltaX
    #print("Derivate              - Backward: " + str(matrix_utils.step_backward_derivate(function, arg1, arg2)))

    #função, x, deltaX
    #print("Derivate              - Central:  " + str(matrix_utils.central_derivate(function, arg1, arg2)))

    #função, x, deltaX, p
    #alc_list7.forward_richard_extrapolation_derivate(function, arg1, arg2, arg3)

    #função, x, deltaX, p
    #alc_list7.backward_richard_extrapolation_derivate(function, arg1, arg2, arg3)

    #função, x, deltaX, p
    #alc_list7.central_richard_extrapolation_derivate(function, arg1, arg2, arg3)


if __name__ == "__main__":

    # Função simples, recebe um único parâmetro x que pode ser uma variável ou um vetor de variáveis
    def function(x):
        return math.log(math.cosh(x*math.sqrt(9.806*0.0034)), 10) - 50

    def function_1(x):
        return 8*x[1]**3 + 6*x[1]*x[0]**2 + 36*x[1]*x[0]*x[2]+108*x[1]*x[2]**2 - 0

    def function_2(x):
        return 60*(x[1]**4) + 60*(x[1]**2)*(x[0]**2)+576*(x[1]**2)*x[0]*x[2]+2232*(x[1]**2)*x[2]**2+252*(x[2]**2)*x[0]**2+1296*(x[2]**2)*x[0]+3348*(x[2]**4)*24*x[0]**3*(x[2])+3*x[0]-3

    # Para cada par de pontos, uma função desse tipo deve ser criada
    def reg_lin_function(x):
        point = [1, 1]
        return x[0] + x[1]*point[0]**x[2] - point[1]

    def reg_lin_function_1(x):
        point = [2, 2]
        return x[0] + x[1]*point[0]**x[2] - point[1]

    def reg_lin_function_2(x):
        point = [3, 9]
        return x[0] + x[1]*point[0]**x[2] - point[1]

    # Este tipo de função tem 2 parâmetros t e y(t)
    def diferential_function(t, y):
        return -2*t*(y**2)

    # Este tipo função tem 3 parâmetros t, y(t) e y'(t)
    def second_order_diferential_function(t, y, z):
        return -9.81 - z*math.sqrt(z**2)

    #main(function, 10, 200, 500)
    #main([function,function_1, function_2], [0.4784, 0.5966, -0.8031])
    main([reg_lin_function, reg_lin_function_1, reg_lin_function_2], [0, 1, 2])
    #main(second_order_diferential_function, 0, 20, 0.1 , 0, 0, False)
