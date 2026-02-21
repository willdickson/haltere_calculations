import sympy
sympy.init_printing()

m = sympy.Symbol('m')
s = sympy.Symbol('s')
t = sympy.Symbol('t')
T = sympy.Symbol('T')
A = sympy.Symbol('A')
Phi = sympy.Symbol('Phi')
f = sympy.Function('f')

if 0:
    expr = A*f.diff(f(t),t)
    series = sympy.fourier_series(expr, (t, -T/2, T/2))
    
    print('expr = ')
    sympy.pprint(expr)
    print()
    
    print('fourier_series = ')
    sympy.pprint(series.truncate(3))
    print()


if 0:
    k = 2
    for n in [1, 3, 5]:
        eq = sympy.sin(2*sympy.pi*t*k/T)*sympy.sin(2*sympy.pi*t/T)*(sympy.cos(2*sympy.pi*t/T)**n)
        #eq = sympy.cos(2*sympy.pi*t*k/T)*sympy.sin(2*sympy.pi*t/T)*(sympy.cos(2*sympy.pi*t/T)**n)
        sol = sympy.integrate(eq, (t, -T/2, T/2))
        print(f'n = {n}')
        sympy.pprint(eq)
        print(sol)
        print()

if 0:
    k = 2
    for n in [1, 3, 5, 7]:
        eq = sympy.sin(2*s)*sympy.sin(s)*sympy.cos(s)**n
        sol = sympy.integrate(eq, (s, -sympy.pi, sympy.pi))
        print(f'n = {n}')
        sympy.pprint(eq)
        print(sol)
        print('-'*30)
        print()




if 0:
    k = 2
    for n in [1, 3, 5, 7]:
        #eq = sympy.cos(s)*sympy.cos(s)**n
        eq = sympy.cos(3*s)*sympy.cos(s)**n
        sol = sympy.integrate(eq, (s, -sympy.pi, sympy.pi))
        print(f'n = {n}')
        sympy.pprint(eq)
        print(sol)
        print('-'*30)
        print()

if 1:
    eq = sympy.cos(s)*sympy.cos(s)**m
    sympy.pprint(eq)
    #eq = sympy.cos(3*s)*sympy.cos(s)**n
    sol = sympy.integrate(eq, (s, -sympy.pi, sympy.pi))
    print(sol)
    print('-'*30)
    print()


