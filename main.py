from importlib.machinery import SourceFileLoader
import math

Matrix_Utils = SourceFileLoader("Matrix_Utils", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/Utils/Matrix_Utils.py").load_module()
ALC_List1    = SourceFileLoader("ALC_List1", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 1/ALC_List1.py").load_module()
ALC_List2    = SourceFileLoader("ALC_List2", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 2/ALC_List2.py").load_module()
ALC_List3    = SourceFileLoader("ALC_List3", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 3/ALC_List3.py").load_module()
ALC_List4    = SourceFileLoader("ALC_List4", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 4/ALC_List4.py").load_module()
ALC_List5    = SourceFileLoader("ALC_List5", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 5/ALC_List5.py").load_module()
ALC_List6    = SourceFileLoader("ALC_List6", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 6/ALC_List6.py").load_module()
ALC_List7    = SourceFileLoader("ALC_List7", "/home/bdantas/Documentos/Repos/Algebra-Linear-Computacional/List 7/ALC_List7.py").load_module()


def main(function, arg1, arg2=[], arg3=[], arg4=[], arg5=[], arg6=[]):
    #Função, menor valor, maior valor
    #ALC_List4.bissection_method(function, arg1, arg2)

    #função, número inicial
    #ALC_List4.newton_method(function, arg1)

    #função, número inicial
    #ALC_List4.secant_method(function, arg1)

    #função, 3 números crescentes próximos à raiz
    #ALC_List4.inverse_interpolation_method(function, arg1, arg2, arg3)

    #lista de funções, chute inicial de solução
    #ALC_List4.non_linear_newton_method(function, arg1)

    #lista de funções, chute inicial de solução
    #ALC_List4.broyden_method(function, arg1)

    #lista de funções reglin, vetor x, vetor y, chute inicial de solução
    ALC_List4.non_linear_mmq(function,arg1)

    #função, a, b
    #ALC_List5.estimate_integral_value(function, arg1, arg2)

    #função, a, b, sensibilidade, boleano (True = polinomial, False = gauss)
    #ALC_List5.integrate(function, arg1, arg2, arg3, False)

    #função diferencial, t0, tf, delta = 0.1, condição inicial, escolha (0 = euler, 1 = RK2, Default = RK4)
    #ALC_List6.first_order_edo_solver(function, arg1, arg2, arg3, arg4, arg5)

    #função diferencial de segunda ordem, t0, tf, delta = 0.1, condição inicial, condição inicial da derivada, escolha(true=Taylor, false=RKN)
    #ALC_List6.second_order_edo_solver(function, arg1, arg2, arg3, arg4, arg5, arg6)

    #função, x, deltaX
    #print("Derivate              - Forward:  " + str(Matrix_Utils.step_forward_derivate(function, arg1, arg2)))

    #função, x, deltaX
    #print("Derivate              - Backward: " + str(Matrix_Utils.step_backward_derivate(function, arg1, arg2)))

    #função, x, deltaX
    #print("Derivate              - Central:  " + str(Matrix_Utils.central_derivate(function, arg1, arg2)))

    #função, x, deltaX, p
    #ALC_List7.forward_richard_extrapolation_derivate(function, arg1, arg2, arg3)

    #função, x, deltaX, p
    #ALC_List7.backward_richard_extrapolation_derivate(function, arg1, arg2, arg3)

    #função, x, deltaX, p
    #ALC_List7.central_richard_extrapolation_derivate(function, arg1, arg2, arg3)

if __name__ == "__main__":

    #Função simples, recebe um único parâmetro x que pode ser uma variável ou um vetor de variáveis
    def function(x):
        return math.log(math.cosh(x*math.sqrt(9.806*0.0034)), 10) - 50

    def function_1(x):
        return 8*x[1]**3 + 6*x[1]*x[0]**2 + 36*x[1]*x[0]*x[2]+108*x[1]*x[2]**2 - 0;

    def function_2(x):
        return 60*(x[1]**4) + 60*(x[1]**2)*(x[0]**2)+576*(x[1]**2)*x[0]*x[2]+2232*(x[1]**2)*x[2]**2+252*(x[2]**2)*x[0]**2+1296*(x[2]**2)*x[0]+3348*(x[2]**4)*24*x[0]**3*(x[2])+3*x[0]-3

    #Para cada par de pontos, uma função desse tipo deve ser criada
    def reg_lin_function(x):
        point = [1, 1]
        return x[0] + x[1]*point[0]**x[2] - point[1]
    def reg_lin_function_1(x):
        point = [2, 2]
        return x[0] + x[1]*point[0]**x[2] - point[1]
    def reg_lin_function_2(x):
        point = [3, 9]
        return x[0] + x[1]*point[0]**x[2] - point[1]

    #Este tipo de função tem 2 parâmetros t e y(t)
    def diferential_function(t, y):
        return -2*t*(y**2)

    #Este tipo função tem 3 parâmetros t, y(t) e y'(t)
    def second_order_diferential_function(t, y, z):
        return -9.81 - z*math.sqrt(z**2)

    #main(function, 10, 200, 500)
    #main([function,function_1, function_2], [0.4784, 0.5966, -0.8031])
    main([reg_lin_function, reg_lin_function_1, reg_lin_function_2], [0,1,2])
    #main(second_order_diferential_function, 0, 20, 0.1 , 0, 0, False)
