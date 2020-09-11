import ALC_List2
import Matrix_Utils

def main(matrix_a, matrix_b=[]):
    print(ALC_List2.power_method(matrix_a, matrix_b))

if __name__ == "__main__":
    main([[3, 2, 0], [2, 3, -1], [0, -1, 3]], [1, -1, 1])
