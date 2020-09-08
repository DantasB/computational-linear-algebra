import ALC_List1

def main(matrix_a, matrix_b=[]):
    print(ALC_List1.solve_by_cholesky_decomposition(matrix_a, matrix_b))

if __name__ == "__main__":
    main([[1, 0.2, 0.4],[0.2, 1, 0.5],[0.4, 0.5, 1]], [0.6, -0.3, -0.6])
