from importlib.machinery import SourceFileLoader

Matrix_Utils = SourceFileLoader("Matrix_Utils", "/home/bdantas/Área de Trabalho/ALC_Lists/Utils/Matrix_Utils.py").load_module()
ALC_List1    = SourceFileLoader("ALC_List1", "/home/bdantas/Área de Trabalho/ALC_Lists/List_1/ALC_List1.py").load_module()
ALC_List2    = SourceFileLoader("ALC_List2", "/home/bdantas/Área de Trabalho/ALC_Lists/List_2/ALC_List2.py").load_module()
ALC_List3    = SourceFileLoader("ALC_List3", "/home/bdantas/Área de Trabalho/ALC_Lists/List_3/ALC_List3.py").load_module()


def main(arg1, arg2=[]):
    print(Matrix_Utils.matrix_determinant(arg1));

if __name__ == "__main__":
    main([[3,1,0],[0,3,0],[0,0,3]], [1,-1,1])
