// Regla de trapecio
function I = trapecio(f, x0, x1)
    I = (x1-x0)*(f(x0)+f(x1))/2
endfunction

function I = trapecioCompuesto(f, a, b, n)
    h = (b-a)/n
    acum = f(a)/2 + f(b)/2
    for j=1:n-1
        acum = acum + f(a+j*h)
    end
    I = h * acum
endfunction

// Regla de Simpson
function I = simpson(f, a, b)
    h = (b-a)/2
    I = h*(f(a)+4*f(a+h)+f(b))/3
endfunction

function I = simpsonCompuesto(f, a, b, n)
    if modulo(n, 2) <> 0 then
        error("simpsonCompuesto - n debe ser un numero par")
        abort()
    end
    h = (b-a)/n
    acum = f(a) + f(b)
    // Los impares
    for j=1:2:n-1
        acum = acum + 4 * f(a + j*h)
    end

    // Los pares
    for j=2:2:n-1
        acum = acum + 2 * f(a + j*h)
    end
    I = (h/3) * acum
endfunction

// Integración numerica en dominio bidimensional

// Funcion que aplica la regla de trapecio para resolver una doble integral
// definida.
function I = trapecioExtendida(f, a, b, c, d)
    // Calculamos h.
    h = (b-a)*(d-c)/4
    // Aplicamos la regla.
    I = h* (f(c,a)+f(c,b)+f(d,a)+f(d,b)) 
endfunction

// Trapecio para f de 2 variables con 'x' fijo
function I = trapecioCompDoblex(f, xi, a, b, n)
    h = (b-a)/n
    contador = f(xi,a)*1/2 + f(xi,b)*1/2
    for j = 1:n-1
        contador = contador + f(xi,a+j*h) 
    end
    I = h*contador
endfunction

// Trapecio para f de 2 variables con 'y' fijo
function I = trapecioCompDobley(f, yi, a, b, n)
    h = (b-a)/n
    contador = f(a,yi)*1/2 + f(b,yi)*1/2
    for j = 1:n-1
        contador = contador + f(a+j*h,yi) 
    end
    I = h*contador
endfunction

// Metdodo de simpson para f de dos variables con 'x' fija
function I = simpsonCompDoblex(f,xi,a,b,n)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        I = %nan
        return
    end
    h = (b-a)/n
    I = f(xi,a)+f(xi,b)
    for i = 1:2:n-1
        I = I+4*f(xi,a+i*h)
    end
    for i = 2:2:n-1
        I = I+2*f(xi,a+i*h)
    end
    I = I*h/3
endfunction

// Metdodo de simpson para f de dos variables con 'y' fija
function I = simpsonCompDobley(f,yi,a,b,n)
    if modulo(n,2) <> 0 then
        disp("N no es par")
        I = %nan
        return
    end
    h = (b-a)/n
    I = f(a,yi)+f(b,yi)
    for i = 1:2:n-1
        I = I+4*f(a+i*h,yi)
    end
    for i = 2:2:n-1
        I = I+2*f(a+i*h,yi)
    end
    I = I*h/3
endfunction

// Metodo del Trapecio para f de 2 variables
// Con dxdy (Primero se integra 'x' luego 'y')
function I = trapecioBidimensionalXY(f, a, b, c, d, n, m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(y) d(y)
    hy = (b-a)/n
    I = trapecioCompDobley(f,a,c(a),d(a),m)/2 + trapecioCompDobley(f,b,c(b),d(b),m)/2
    for i = 1:n-1
        yi = a+i*hy
        I = I + trapecioCompDobley(f,yi,c(yi),d(yi),m)
    end
    I = hy*I
endfunction

// Metodo del Trapecio para f de 2 variables
// Con dydx (Primero se integra 'y' luego 'x')
function trapecioBidimensionalYX(f, a, b, c, d, n, m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(x) d(x)
    hx = (b-a)/n
    I = trapecioCompDoblex(f,a,c(a),d(a),m)/2 + trapecioCompDoblex(f,b,c(b),d(b),m)/2
    for i = 1:n-1
        xi = a+i*hx
        I = I + trapecioCompDoblex(f,xi,c(xi),d(xi),m)
    end
    I = hx*I
endfunction

// Metodo de Simpson para f de 2 variables
// Con dxdy (Primero se integra 'x' luego 'y')
function I = simpsonCompBidimensionalXY(f, a, b, c, d, n, m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(y) d(y)
    // m y n deben ser par
    if modulo(n,2) <> 0 then
        disp("N no es par")
        I = %nan
        return
    end
    hy = (b-a)/n
    I = simpsonCompDobley(f,a,c(a),d(a),m) + simpsonCompDobley(f,b,c(b),d(b),m)
    for i = 1:2:n-1
        yi = a+i*hy
        I = I + 4*simpsonCompDoblex(f,yi,c(yi),d(yi),m)
    end
    for i = 2:2:n-1
        yi = a+i*hy
        I = I + 2*simpsonCompDoblex(f,yi,c(yi),d(yi),m)
    end
    I = I*hy/3
endfunction

// Metodo de Simpson para f de 2 variables
// Con dydx (Primero se integra 'y' luego 'x')
function I = simpsonCompBidimensionalYX(f, a, b, c, d, n, m)
    // a,b valores de integral exterior
    // c,d valores de integral interior (Se reciben como funcion, si se quiere constante ingresar funciones constantes)
    // n intervalos entre a y b
    // m intervalos en c(x) d(x)
    // m y n deben ser par
    if modulo(n,2) <> 0 then
        disp("N no es par")
        I = %nan
        return
    end
    hx = (b-a)/n
    I = simpsonCompDoblex(f,a,c(a),d(a),m) + simpsonCompDoblex(f,b,c(b),d(b),m)
    for i = 1:2:n-1
        xi = a+i*hx
        I = I + 4*simpsonCompDoblex(f,xi,c(xi),d(xi),m)
    end
    for i = 2:2:n-1
        xi = a+i*hx
        I = I + 2*simpsonCompDoblex(f,xi,c(xi),d(xi),m)
    end
    I = I*hx/3
endfunction

// METODOS INTEGRACION BIDIMENSIONAL DADOS EN CLASE

// Funcion que toma una funcion, un intervalo y una cantidad de subIntervalos.
function res = metodoCompTrapecio(fx, a, b, n)
    // Inicializamos h y res
    h = (b-a)/n;
    res = 0;
    // Iteramos sobre todos los subintervalos.
    for i = 0:n
        xn = a + i*h;
        // Si el intervalo es el primero o el ultimo dividimos f(xn) por 2.
        if(i == 0 | i == n)
            res = res + fx(xn)/2;
        // En caso contrario no.
        else
            res = res + fx(xn);
        end
    end
    // Multiplicamos todo por h.
    res = res * h;
endfunction


// Fucnion que calcula el metodo compuesto de Simpson.
function res = metodoCompSimpson(fx, a, b, n)
    // Inicializamos h y res.
    h = (b-a)/n
    res = 0
    
   // Primer intervalo.
    res = res + fx(a + 0*h)
    // Iteramos sobre los n subintervalos.
    for i = 1:n-1
        
        // Si el subintervalo es par.
        if (pmodulo(i,2) == 0)
            res = res + 2*fx(a + i*h);
        // Caso impar.
        else
            res = res + 4*fx(a + i*h);
        end    
        
    end
    // Ultimos intervalo.
    res = res + fx(a + n*h)
    // Multiplicamos por h/3
    res = res * h/3;
endfunction

// Funcion que calcula la integral numerica de dominio bidimensional, 
// utilizando el metodo de Simpson.
// Toma la funcion, los extremos del primer intervalo, las funciones del segundo
// y la cantidad de subintervalos de la primer integral y la cantidad de subintervalos
// de la segunda.
function y = dominioBiSimpson(f,a,b,cx,dx,n,m)
    // Inicializamos h
    h = (b-a)/n
    deff('z=fxa(y)','z=f('+string(a)+',y)')
    deff('z=fxb(y)','z=f('+string(b)+',y)')
    
    
    // Calculamos la segunda integral.
    temp = metodoCompSimpson(fxa,cx(a),dx(a),m) + metodoCompSimpson(fxb,cx(b),dx(b),m)
    
    // Calculamos la primera, en todos los subintervalos pedidos.
    for i=1:n-1
        xi = a+i*h
        deff('z=aux(y)','z=f(xi,y)')
        // Aplicamos el metodo como corresponde segun simpson.
        if pmodulo(i,2) == 0 then
            temp = temp + 2*(metodoCompSimpson(aux,cx(xi),dx(xi),m))
        else
            temp = temp + 4*(metodoCompSimpson(aux,cx(xi),dx(xi),m))
        end
    end
    y = (h/3) * temp
endfunction

// Funcion que calcula la integral numerica de dominio bidimensional, 
// utilizando el metodo de Trapecio.
// Toma la funcion, los extremos del primer intervalo, las fucniones del segundo
// y la cntidad de subintervalos de la primer integral y la cantidad de subintervalos
// de la segunda.
function y = dominioBiTrapecio(f,a,b,cx,dx,n,m)
    // Inicializamos h
    h = (b-a)/n
    
    deff('z=fxa(y)','z=f('+string(a)+',y)')
    deff('z=fxb(y)','z=f('+string(b)+',y)')
    
    
    // Calculamos la segunda integral.
    temp = (metodoCompTrapecio(fxa, cx(a),dx(a),m)/2) + (metodoCompTrapecio(fxb,cx(b),dx(b),m)/2)
    
    // Calculamos la primera, en todos los subintervalos pedidos.

    for i=1:n-1
        xi = a+i*h
        deff('z=aux(y)','z=f('+string(xi)+',y)')
        temp = temp + (metodoCompTrapecio(aux,cx(xi),dx(xi),m))
    end
    y = h * temp
endfunction

// Ejercicios

// Ejercicio 1
// log = ln

// apartado a
/*
inf = 1
sup = 2
I1 = trapecio(log, inf, sup)
I2 = simpson(log, inf, sup)
valorReal = intg(inf, sup, log)

disp("El valor por trapecio es: ")
disp(I1)
disp("El valor por simpson es: ")
disp(I2)
disp("El valor real segun scilab es: ")
disp(valorReal)

// apartado b

inf = 0
sup = 0.1

function y = fb(x)
    y = x^(1/3)
endfunction

I1 = trapecio(fb, inf, sup)
I2 = simpson(fb, inf, sup)
valorReal = intg(inf, sup, fb)

disp("El valor por trapecio es: ")
disp(I1)
disp("El valor por simpson es: ")
disp(I2)
disp("El valor real segun scilab es: ")
disp(valorReal)

// apartado c

function y = fc(x)
    y = sin(x)^2
endfunction

inf = 0
sup = %pi/3
I1 = trapecio(fc, inf, sup)
I2 = simpson(fc, inf, sup)
valorReal = intg(inf, sup, fc)

disp("El valor por trapecio es: ")
disp(I1)
disp("El valor por simpson es: ")
disp(I2)
disp("El valor real segun scilab es: ")
disp(valorReal)
*/
// apartado iii en papel

// Ejercicio 2 y 3
/*
// a

function y = f2a(x)
    y = 1/x
endfunction

a = 1
b = 3
n = 4
ITrap = metodoCompTrapecio(f2a, a, b, n)
ISimp = metodoCompSimpson(f2a, a, b, n)

valor = intg(a, b, f2a)

disp("Apartado a")

printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))

// b

function y = f2b(x)
    y = x^3
endfunction

a = 0
b = 2
n = 4

ITrap = metodoCompTrapecio(f2b, a, b, n)
ISimp = metodoCompSimpson(f2b, a, b, n)

valor = intg(a, b, f2b)

disp("Apartado b")

printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))

// c

function y = f2c(x)
    y = x * ((1+x^2)^1/2)
endfunction

a = 0
b = 3
n = 6
ITrap = metodoCompTrapecio(f2c, a, b, n)
ISimp = metodoCompSimpson(f2c, a, b, n)

valor = intg(a, b, f2c)

disp("Apartado c")
printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))

// d

function y = f2d(x)
    y = sin(%pi*x)
endfunction

a = 0
b = 1
n = 8
ITrap = metodoCompTrapecio(f2d, a, b, n)
ISimp = metodoCompSimpson(f2d, a, b, n)

valor = intg(a, b, f2d)

disp("Apartado d")
printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))

// e

function y = f2e(x)
    y = x * sin(x)
endfunction

a = 0
b = 2*%pi
n = 8
ITrap = metodoCompTrapecio(f2e, a, b, n)
ISimp = metodoCompSimpson(f2e, a, b, n)

valor = intg(a, b, f2e)

disp("Apartado e")

printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))

// f

function y = f2f(x)
    y = (x^2)*(%e^x)
endfunction

a = 0
b = 1
n = 8
ITrap = metodoCompTrapecio(f2f, a, b, n)
ISimp = metodoCompSimpson(f2f, a, b, n)

valor = intg(a, b, f2f)

disp("Apartado f")

printf("Valor exacto = %f\n", valor)
disp("Metodo compuesto Trapecio")
disp(ITrap)
disp("Error Trapecio")
disp(abs(ITrap - valor))
disp("Metodo compuesto Simpson")
disp(ISimp)
disp("Error Simpson")
disp(abs(ISimp - valor))
*/
// Ejercicio 4
/*
function y = f(x)
    y = (x + 1)^-1
endfunction

a = 0
b = 1.5
n = 10
valor = 0.9262907
I1 = metodoCompTrapecio(f, a, b, n)
I2 = metodoCompSimpson(f, a, b, n)
disp("Valor exacto = 0.9262907")

disp("Metodo compuesto Trapecio")
disp(I1)
disp("Error Trapecio")
disp(abs(I1 - valor))
disp("Metodo compuesto Simpson")
disp(I2)
disp("Error Simpson")
disp(abs(I2 - valor))
*/

// Ejercicio 5
/*
function z=f5(x,y)
    z = sin(x+y)
endfunction

disp("La aproximacion de la integral doble por trapecio extendida es")
disp(trapecioExtendida(f5,0,2,0,1))
*/
// Ejercicio 6
/*
function y = cosenoinf1(x)
    y = -sqrt(2*x - x**2)
endfunction
function y = cosenosup1(x)
    y = sqrt(2*x - x**2)
endfunction

function y = cosenoinf2(x)
    y = -sqrt(1- x**2)
endfunction
function y = cosenosup2(x)
    y = sqrt(1 - x**2)
endfunction

function z = circulo(x,y)
    z = 1
endfunction

inf = 0
sup = 2
n = 2 // 150 ya tiene una buena aproximacion de pi
m = 2 // 150 ya tiene una buena aproximacion de pi

I1 = dominioBiTrapecio(circulo,inf,sup,cosenoinf1,cosenosup1,n,m)
I2 = dominioBiTrapecio(circulo,-1,1,cosenoinf2,cosenosup2,n,m)
I3 = dominioBiSimpson(circulo,inf,sup,cosenoinf1,cosenosup1,n,m)
I4 = dominioBiSimpson(circulo,-1,1,cosenoinf2,cosenosup2,n,m)

printf("Aproximacion area circulo por trapecio compuesto bidimensional YX: %.10f\n", I1)
printf("Aproximacion area circulo por trapecio compuesto bidimensional XY: %.10f\n", I2)
printf("Aproximacion area circulo por simpson compuesto bidimensional YX: %.10f\n", I3)
printf("Aproximacion area circulo por simpson compuesto bidimensional XY: %.10f\n", I4)
printf("Pi == %.10lf\n",%pi)
printf("error I1: %.10f\n",abs(I1-%pi))
printf("error I2: %.10f\n",abs(I2-%pi))
printf("error I3: %.10f\n",abs(I3-%pi))
printf("error I4: %.10f\n",abs(I4-%pi))

// Si cambiamos los n y m podemos observar como al aumentar los subintervalos el resultado se acerca a pi

*/

// Ej 6 v2
/*
function y=f(x,y)
    y=1
endfunction

function y=dx(x)
    y= sqrt(2*x-x^2)
endfunction

function y=cx(x)
    y= -sqrt(2*x-x^2)
endfunction

disp(dominioBiTrapecio(f,0,2,cx,dx,2,2))
// 2

disp(dominioBiTrapecio(f,0,2,cx,dx,100,100))
// 3.1382685

disp(dominioBiTrapecio(f,0,2,cx,dx,1000,1000))
// 3.1414875


disp(dominioBiSimpson(f,0,2,cx,dx,2,2))
// 2.6666667

disp(dominioBiSimpson(f,0,2,cx,dx,100,100))
// 3.140296

disp(dominioBiSimpson(f,0,2,cx,dx,1000,1000))
// 3.1415516

// Vemos como al aumentar los subintervalos el resultado se acerca a pi
*/
