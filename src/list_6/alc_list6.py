import matplotlib.pyplot as plt
import math


def first_order_edo_solver(diferential_function, t0, tf, delta, start_condition, control):
    x_incognita = [start_condition]
    t_incognita = [t0]
    steps = int((tf - t0) / (delta))

    if(control == 0):
        print("Euler Method Solution")
    elif(control == 1):
        print("Runge Kutta Second Order Method Solution")
    else:
        print("Runge Kutta Fourth Order Method Solution")

    print(" t", "     ", " x")
    print("0.0", "    ", start_condition)

    for i in range(steps):
        t_incognita.append((i + 1) * delta)

        if(control == 0):
            x_incognita.append(
                x_incognita[i] + delta * diferential_function(t_incognita[i], x_incognita[i]))
        elif(control == 1):
            K1 = diferential_function(t_incognita[i], x_incognita[i])
            K2 = diferential_function(
                t_incognita[i] + delta, x_incognita[i] + delta * K1)
            x_incognita.append(x_incognita[i] + delta / 2 * (K1 + K2))
        else:
            K1 = diferential_function(t_incognita[i], x_incognita[i])
            K2 = diferential_function(
                t_incognita[i] + delta / 2, x_incognita[i] + delta / 2 * K1)
            K3 = diferential_function(
                t_incognita[i] + delta / 2, x_incognita[i] + delta / 2 * K2)
            K4 = diferential_function(
                t_incognita[i] + delta, x_incognita[i] + delta * K3)
            x_incognita.append(
                x_incognita[i] + delta / 6 * (K1 + 2 * K2 + 2 * K3 + K4))

        if(i == 0 or i == steps):
            continue

        print(round(t_incognita[i], 3), "    ", round(x_incognita[i], 5))

    print(round(t_incognita[-1], 3), "   ", round(x_incognita[-1], 5))
    print('\n')


def second_order_edo_solver(second_order_diferential_function, t0, tf, delta, x0, derivate_x0, control=True):
    x_incognita = [x0]
    x_actual_line = derivate_x0
    t_incognita = [t0]
    steps = int((tf - t0) / delta)
    for i in range(steps):
        t_incognita.append((i + 1) * delta)

        if(control):
            x2_lines = second_order_diferential_function(
                t_incognita[i], x_incognita[i], x_actual_line)
            x_incognita.append(
                x_incognita[i] + x_actual_line * delta + x2_lines * delta * delta / 2)
            x_actual_line = x_actual_line + x2_lines * delta

        else:
            K1 = delta / 2 * \
                second_order_diferential_function(
                    t_incognita[i], x_incognita[i], x_actual_line)
            Q = delta / 2 * (x_actual_line + K1 / 2)
            K2 = delta / 2 * second_order_diferential_function(
                t_incognita[i] + delta / 2, x_incognita[i] + Q, x_actual_line + K1)
            K3 = delta / 2 * second_order_diferential_function(
                t_incognita[i] + delta / 2, x_incognita[i] + Q, x_actual_line + K2)
            L = delta * (x_actual_line + K3)
            K4 = delta / 2 * second_order_diferential_function(
                t_incognita[i] + delta, x_incognita[i] + L, x_actual_line + 2*K3)

            x_incognita.append(
                x_incognita[i] + delta*(x_actual_line + (K1 + K2 + K3) / 3))
            x_actual_line = x_actual_line + (K1 + 2 * K2 + 2 * K3 + K4) / 3

    plt.plot(t_incognita, x_incognita)
    plt.ylabel('y\'\'(t)')
    plt.xlabel('T')
    plt.show()
