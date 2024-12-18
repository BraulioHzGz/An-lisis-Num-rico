'''
    CÓDIGOS PARA LA MATERIA DE ANÁLISIS NUMÉRICO
    4o SEMESTRE DE LA FACULTAD DE INGENIERÍA UNAM

    CÓDIGOS FUENTES: PYTHON (.py)
    REQUERIMENTOS: INSTALAR PAQUETES NECESARIOS
'''

import numpy as np
import sympy as sym
import math
import warnings
from tabulate import tabulate
warnings.filterwarnings("ignore")

def ErrorAbsoluto(exacto, aproximado):
    return abs(exacto - aproximado)


def ErrorRelativo(exacto, aproximado):
    return (ErrorAbsoluto(exacto, aproximado) / abs(exacto))


def ErrorRelativoPorcentual(exacto, aproximado):
    return ErrorRelativo(exacto, aproximado) * 100


def PolinomioTaylor(punto, terminos, redondeo):
    print('\nMÉTODO DEL POLINOMIO DE TAYLOR\n')
    # Forzar a 2 decimales de preferencia | Invocarlo como PolinomioTaylor(x0, k, redondeo=2)

    x = sym.Symbol('x')
    fx = sym.cos(x)     # Aquí le cambian f(x) a su necesidad
    x0 = punto
    grado = terminos
    k = 0
    polinomio = 0
    while(k <= grado):
        dfx = fx.diff(x,k)
        dfx0 = dfx.subs(x,x0)
        divisor = math.factorial(k)
        tK = round((dfx0 / divisor), redondeo) * (x - x0)**k
        polinomio = polinomio + tK
        k += 1
    
    print(polinomio)


# Aquí escriben el polinomio que quieran probar en los distintos métodos numéricos
def f(x):
    return x**3 + 3*(x**2) -1


def MetodoBiseccion(a, b):
    print('\nMÉTODO DE LA BISECCIÓN\n')
    print("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}".format('It', 'a', 'b', 'c', 'f(a)', 'f(b)', 'f(c)', 'Ea'))
    it = 1
    toleracia = 0.001   # Entre más pequeño mejor
    while True:
        c = (a+b) / 2
        fa = f(a)
        fb = f(b)
        fc = f(c)
        Ea = ErrorAbsoluto(a, b)
        print("{:<8} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g}".format(it, a, b, c, fa, fb, fc, Ea))
        it += 1
        if(Ea < toleracia):
            print(f'\nRaíz encontrada: {c}')
            break
        elif(fa * fc) < 0: b = c
        elif(fb * fc) < 0: a = c
        else:
            print(f'\nNo hay raíz encontrada en el intervalo que va de {a} a {b}')
            break


def ReglaFalsa(a, b):
    print('\nMÉTODO DE LA REGLA FALSA\n')
    tolerancia = 0.001
    print("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} ".format('It', 'a', 'b', 'c', 'f(a)', 'f(b)', 'f(c)'))
    it = 1
    while True:
        it += 1
        c = b - ((f(b) * (a - b)) / (f(a) - f(b)))
        fa = f(a)
        fb = f(b)
        fc = f(c)
        print("{:<8} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g} {:<8.4g}".format(it, a, b, c, fa, fb, fc))
        if(abs(fc) < tolerancia):
            print(f'\nRaíz encontrada: {c}')
            break
        elif(fa * fc) < 0: b = c
        elif(fb * fc) < 0: a = c
        else:
            print(f'\nNo hay raíz encontrada en el intervalo que va de {a} a {b}')
            break


def NewtonRhapson(ValorInicial):
    print('\nMÉTODO DE NEWTON - RHAPSON\n')
    x = sym.Symbol('x')
    f_sym = x**2 - 10 * sym.cos(x)
    df_sym = f_sym.diff(x)

    fx = sym.lambdify(x, f_sym, modules='math')
    dfx = sym.lambdify(x, df_sym, modules='math')
    xi = ValorInicial
    toleracia = 0.001
    it = 1
    print("{:<8} {:<8} {:<8} {:<8} ".format('It', 'x_i', 'x_i+1', 'Ea'))
    while True:
        if(dfx(xi) == 0):
            print("\nNo se encontró una raíz")
            break
        xi_prime = xi - fx(xi) / dfx(xi)
        Ea = ErrorAbsoluto(xi_prime, xi)
        print("{:<8} {:<8.4g} {:<8.4g} {:<8.4g} ".format(it, xi, xi_prime, Ea))
        xi = xi_prime
        it += 1
        if(Ea < toleracia):
            print(f'\nRaíz encontrada: {xi_prime}')
            break


def MetodoLin(coeficientes, decimales):
    print('\nMÉTODO DE LIN\n')
    grado = len(coeficientes) - 1
    if grado < 3:
        print("\nEl polinomio debe ser de grado 3 o mayor")
        return
    
    p = q = 0
    tolerancia = 0.001

    NuevosValores = [p, q, 0, 0]    # [p, q, Δp, Δq]
    if grado > 3: table = [["p", "q", "b0", "b1", "b2", "R", "S", "Δp", "Δq"]]
    else: table = [["p", "q", "b0", "b1", "R", "S", "Δp", "Δq"]]
    
    def CalcularValores(grado, valores):
        p, q, IncP, IncQ = valores
        b0 = coeficientes[0]
        b1 = round(coeficientes[1] - (p * b0), decimales)
        b2 = round(coeficientes[2] - (p * b1) - (q * b0), decimales) if grado > 3 else None
        R = round(coeficientes[3] - (p * b2) - (q * b1), decimales) if grado > 3 else round(coeficientes[2] - (p * b1) - (q * b0), decimales)
        S = round(coeficientes[4] - (q * b2), decimales) if grado > 3 else round(coeficientes[3] - (q * b1), decimales)
        IncP = round(R / (b2 if grado > 3 else b1), decimales)
        IncQ = round(S / (b2 if grado > 3 else b1), decimales)

        fila = [p, q, b0, b1, b2, R, S, IncP, IncQ] if grado > 3 else [p, q, b0, b1, R, S, IncP, IncQ]
        return [p + IncP, q + IncQ, IncP, IncQ], fila

    NuevosValores, fila = CalcularValores(grado, NuevosValores)
    table.append(fila)

    while(abs(NuevosValores[2]) > tolerancia or abs(NuevosValores[3]) > tolerancia):
        NuevosValores, fila = CalcularValores(grado, NuevosValores)
        table.append(fila)

    print(tabulate(table, tablefmt='fancy_grid', stralign='center'))


def GaussJordan(matrix):
    print('\nMÉTODO DE GAUSS JORDAN\n')
    print("\nSistema a resolver:\n\na1\t\ta2\t\ta3\t\tcte.")
    for i in matrix:
        for j in i: print(j, end='\t\t')
        print("\n")
    
    def Ceros(a, b):
        for i in range(len(matrix[0])):
            if matrix[a][b] != 0:
                c = matrix[a][b]
                for j in range(len(matrix[0])): matrix[a][j] = matrix[a][j] - (c * matrix[b][j])
    
    def Unos(x):
        for i in range (len(matrix[0])):
            if matrix[x][x] != 1:
                y = matrix[x][x]
                for j in range(len(matrix[0])): matrix[x][j] = matrix[x][j] / y
    
    for i in range(len(matrix)):
        Unos(i)
        for j in range(len(matrix)):
            if(i != j): Ceros(j, i)

    print("\nSistema a resuelto:\n\na1\t\ta2\t\ta3\t\tcte.")
    for i in matrix:
        for j in i: print(j, end='\t\t')
        print("\n")


def Gauss_Seidel(A, B):
    print('\nMÉTODO DE GAUSS SEIDEL\n')
    A = np.array(A, dtype = np.float32)
    B = np.array(B, dtype = np.float32)
    it = 1
    tolerancia = 0.1
    n = len(A[0:])
    Z = np.zeros((n))
    print("{:<8} {:<20} {:<8} ".format('It', 'Solución', '      Ea'))

    while True:
        xi = Z[0]
        for i in range(n):
            pivote = A[i, i]
            for j in range(n): A[i, j] = A[i, j] / pivote
            B[i] = B[i] / pivote
        for i in range(n):
            sum = B[i]
            for j in range(n): 
                if(i != j): sum -= - A[i, j] * Z[j]
            Z[i] = sum
        Ea = np.abs(xi - Z[0])
        Z_str = "  ".join([f'{z:.4g}' for z in Z])
        print("{:<8} {:<20} {:<8} ".format(it, "[" + str(Z_str) + "]", "   " + str(np.max(Ea))))
        it += 1
        if(np.max(Ea) < tolerancia): break


def MetodoLu(A, B):
    print('\nMÉTODO DE LU\n')
    A = np.array(A, dtype = np.float32)
    B = np.array(B, dtype = np.float32)
    n = len(B)
    x = y = np.zeros([n, 1])
    U = np.zeros([n, n])
    L = np.identity(n)

    for i in range(n):
        for j in range(i + 1, n):
            factor = (A[j, i] / A[i, i])
            L[j, i] = factor
            for c in range(n): A[j, c] = A[j, c] - factor * A[i, c]
        
    # Obteniendo la solución para y
    y[0] = B[0]
    for i in range(1, n):
        suma = 0
        for j in range(n): suma += (L[i, j] * y[j])
        y[i] = B[i] - suma

    U = A
    x[n - 1] = y[n - 1] / U[n - 1, n - 1]
    for i in range(n - 2, -1, -1):
        suma = 0
        for j in range(n): suma += (U[i, j] * x[j])
        x[i] = (y[i] - suma) / U[i, i]
    
    print('Resultado para L: \n', L)
    print('\nResultado para U: \n', U)
    print('\nResultado para x: \n', x)
    print('\nResultado para y: \n', y)
    print('\nComprobando A*x : \n', A*x)
