function p = lagrange(x,y)
    [nx, mx] = size(x)
    p = 0
    for k=1:mx
        p = p + (Lk(x,k)*y(k))
    end 
endfunction

function l = Lk(x,k)
    [nx, mx] = size(x)
    // Se genera un vector con las raices del polinomio Lk
    r = [x(1:k-1) x(k+1:mx)]
    // Con poly generamos el polinomio cuyas raices son r
    p = poly(r,"x","roots")
    // Evaluamos
    pk = horner(p, x(k))
    l = p / pk
endfunction

// Error de interpolación
function w = err(p,x,cot)
    // Entrada: p = valor real, x = nodos de interpolación, cot = cota de |f^(n))|
    // Salida: w = error de interpolación en x = p
    n = length(x)
    w = cot/(factorial(n+1))
    for i=1:n do
        w = w*abs(p - x(i))
    end
endfunction

// Newton

// una forma de diferencias divididas
function d = difDiv(x, y, i, k)
    if k - i == 1 then
        d = (y(k) - y(i)) / (x(k) - x(i))
    else
        d = (difDiv(x, y, i+1, k) - difDiv(x, y, i, k-1)) / (x(k) - x(i))
    end
endfunction

// otra forma
function d = difDividida(x,y)
    [nx,mx] = size(x)
    if mx == 1 then
        d = y(1)
    else
        d = (difDividida(x(2:mx),y(2:mx))- difDividida(x(1:mx-1),y(1:mx-1)))/(x(mx)-x(1))
    end
endfunction

// con difDiv
function p = newton(x, y)
    // Obtenemos y guardamos el tamaño
    [Xn, Xm] = size(x)
    // calculamos el primer coeficiente
    p = y(1)
    for i=1:Xm-1
        // formula del polinomio de newton.
        p = p + poly(x(1:i),"x","roots") * difDiv(x, y, 1, i+1)
    end
endfunction

//con difDividida
function p = newton2(x,y)
    [nx,mx] = size(x)
    p = y(1)
    for k=2:mx
        q = poly(x(1:k-1),"x","roots")
        dif = difDividida(x(1:k),y(1:k))
        p = p + dif*q
    end  
endfunction

// Diferencias divididas
function w = DD(x,y)
    // Entrada: x,y = vectores de puntos de interpolación (x,y)
    // Salida: w = diferencias divididas en los vectores x e y
    n = length(x)
    if n==2 then
        w = (y(n)-y(1))/(x(n)-x(1))
    else
        w = (DD(x(2:n),y(2:n))-DD(x(1:n-1),y(1:n-1)))/(x(n)-x(1))
    end
endfunction

// Método de Diferencias Divididas de Newton
// Formula con multiplicaciones encajadas
function w = DD_NewtonC(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    s = poly(0,"x")
    n = length(x)
    w =  difDividida(x,y)
    for j=n-1:-1:1
        w = (s-x(j))*w + DD(x(1:j),y(1:j))
    end
endfunction

// Método de Diferencias Divididas de Newton
function w = DD_Newton(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    w = 0
    s = poly(0,"x")
    n = length(x)
    for j=n:-1:2
        w = w*(s-x(j-1)) + DD(x(1:j),y(1:j))*(s-x(j-1))
    end
    w = w + y(1)
endfunction

// Ejercicio 1
x = [0 0.2 0.4 0.6]
y = [1 1.2214 1.4918 1.8221]
resultado = 1.395612425
v = 1/3

// Lagrange

lagrangeLineal = lagrange(x(2:3),y(2:3))
disp("Pol de lagrange lineal")
disp(lagrangeLineal)
disp("Lagrange evaluado en 1/3")
disp(horner(lagrangeLineal,v))
// 1.4016667
disp("Cota de error lineal:")
disp(abs(((1/3)-0.2)*((1/3)-0.4)*%e^(0.6))/2)
// 0.0080983
disp("Error del polinomio lineal")
disp(abs(resultado-horner(lagrangeLineal,1/3)))
//0.0060542

// cubica

lagrangeCubica = lagrange(x,y)
disp("Pol de lagrange cubico")
disp(lagrangeCubica)
disp("Lagrange evaluado en 1/3")
disp(horner(lagrangeCubica,v))
// 1.3955494
disp("Cota de error cubico:")
disp(abs(((1/3)-0)*((1/3)-0.2)*((1/3)-0.4)*((1/3)-0.6)*%e^(0.6))/factorial(4))
// 0.00006
disp("Error del polinomio cubico")
disp(abs(resultado-horner(lagrangeCubica,1/3)))
// 0.000063

// Newton 

newtonLineal = newton(x(2:3),y(2:3))
newtonCubica = newton(x,y)

disp("Newton")

// Lineal

disp("Resultado de evaluar polinomio lineal en 1/3")
disp(horner(newtonLineal,1/3))
// = 1.4016667
disp("Cota de error lineal:")
disp(abs(((1/3)-0.2)*((1/3)-0.4)*%e^(0.6))/2)
// = 0.0080983
disp("Error del polinomio lineal")
disp(abs(1.395612425-horner(newtonLineal,1/3)))
// = 0.0060542

// Cubico

disp("Resultado de evaluar polinomio cubico en 1/3")
disp(horner(newtonCubica,1/3))
// = 1.3955494
disp("Cota de error cubico:")
disp(abs((1/3-0)*(1/3-0.2)*(1/3-0.4)*(1/3-0.6)*%e^0.6)/factorial(4))
// = 0.00006
disp("Error del polinomio cubico")
disp(abs(1.395612425-horner(newtonCubica,1/3)))
// = 0.000063

// Graficas
rango = [-2:0.01:2]
plot(rango, horner(lagrangeLineal, rango),"r")
plot(rango, horner(lagrangeCubica, rango),"g")
plot(rango, exp(rango))
a=gca();a.x_location = "origin";a.y_location = "origin"
h1 = legend(['lagrangeLineal','lagrangeCubico','e^x'])
scf()

rango = [-2:0.01:2]
plot(rango, horner(newtonLineal, rango),"r")
plot(rango, horner(newtonCubica, rango),"g")
plot(rango, exp(rango))
a=gca();a.x_location = "origin";a.y_location = "origin"
h1 = legend(['newtonLineal','newtonCubica','e^x'])

// para maximizar phi paso las raices poly(raices, "x", "roots") despues derivo con derivat(ese polinomio) y busco las raices del derivat roots(derivat(ese polinomio))
// despues simplemente por fuerza bruta horner(poly de phi, raiz_iesima) y me quedo el max abs
// para maximizar la f^n+1 hago f^n+2 busco sus ceros y me fijo en los extremos en esos puntos cuanto da f^n+1 y me quedo con el mayor valor abs

// segun todo el internet se acota el f^n+1 y listo pero idk
// aparentemente maximizar ambos es para general en toda x y tomando el x = punto en el que evaluas y maximizar f^n+1 va

// Ejercicio 4
x = [2.0 2.1 2.2 2.3 2.4 2.5]
y = [0.2239 0.1666 0.1104 0.0555 0.0025 -0.0484]

p = DD_Newton(x,y)
disp("El polinomio interpolante es: ")
disp(p)

w1 = horner(p,2.15)
err1 = err(2.15,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.15) es: "+string(w1))
disp("con error: "+string(err1)+" < 0.5D-06")

w2 = horner(p,2.35)
err2 = err(2.35,x,1) // |j_0'(x)| <= 1
disp("El valor aproximado de J_0(2.35) es: "+string(w2))
disp("con error: "+string(err2)+" < 0.5D-06")

// Ajuste de curvas

// Es necesario dividir el proceso de la matriz y el polinomio resultado
// ya que se debe checkear que la matriz por su transpuesta sea rango completo
function A = matriz_mc(x,n)
    [nx,mx] = size(x)
    A = ones(mx,1)
    for i=1:n
        A = [A (x')^i]
    end
endfunction

// Aproximación polinomial de mínimos cuadrados polinomial para matrices con rango completo
function [p, er] = MinCuad_pol(A,b)
    // Entrada: b = vectores 1xn
    // Salida: p = polinomio de mínimos cuadrados; err = vector de errores (eps = Ax-b)
     [w,a] = gausselimPP((A')*A,(A')*(b'))
     p = poly(w,"x","coeff")
     er = A*w-b'
endfunction

function [p, erV, A, erN] = minimosCuadrados(x, y, n)
    // salida: 
    // p = el polinomio de minimos cuadrados;
    // erV = vector de errores;
    // A = matriz del método de mínimo cuadrados;
    // erN = el error (norm(Ax-y,2))
    A = matriz_mc(x, n)
    [p, erV] = MinCuad_pol(A, y)
    erN = norm(erV, 2)
endfunction

function [x,a] = gausselimPP(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselimPP - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselimPP - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada
n = nA;    // Tamaño de la matriz

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(a(k,k));  //pivoteo
    for i=k+1:n
        if abs(a(i,k))>amax then
            kpivot = i; amax = a(i,k);
        end;
    end;
    temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
    
    for i=k+1:n
        for j=k+1:n+1
            a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
        end;
        for j=1:k        // no hace falta para calcular la solución x
            a(i,j) = 0;  // no hace falta para calcular la solución x
        end              // no hace falta para calcular la solución x
    end;
end;

// Sustitución regresiva
x(n) = a(n,n+1)/a(n,n);
for i = n-1:-1:1
    sumk = 0
    for k=i+1:n
        sumk = sumk + a(i,k)*x(k);
    end;
    x(i) = (a(i,n+1)-sumk)/a(i,i);
end;
endfunction

// Recibe un numero n
// Devuelve el polinomio de chebyshev de ese grado con sus raíces
function [T,r] = Chebyshev(n)
    t(1) = 1
    t(2) = poly([0],"x","r")
    for i = 3:n+1
        t(i) = poly([0 2], "x", "coeff")*t(i-1)-t(i-2) 
    end
    T = t(n+1)
    r = roots(T)
endfunction

// Raices de Chebyshev en cualquier intervalo
function x = NodosChebyshev(n,a,b)
    [pol,r] = Chebyshev(n)
    for i = 1 : n
        x(i) = ((b+a) + r(i) * (b - a))/2
    end
endfunction

// Funcion que calcula los nodos de raices segun el grado.
function r = Cheb(n)
    for k=0:n-1
        r(k+1) = cos(%pi/2*(1+2*k)/n)
    end
endfunction

// Ejercicio 7
x = [0 0.15 0.31 0.5 0.6 0.75]
y = [1 1.004 1.31 1.117 1.223 1.422]
// grado 1
[p1, er] = minimosCuadrados(x, y, 1)
disp("El polinomio de grado 1 es: ")
disp(p1)
// grado 2
[p2, er] = minimosCuadrados(x, y, 2)
disp("El polinomio de grado 2 es: ")
disp(p2)
// grado 3
[p3, er] = minimosCuadrados(x, y, 3)
disp("El polinomio de grado 3 es: ")
disp(p3)
rango = [-0:0.0001:1]
scf()
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])

// Ejercicio 8
x = [4 4.2 4.5 4.7 5.1 5.5 5.9 6.3 6.8 7.1]
y = [102.56 113.18 130.11 142.05 167.53 195.14 224.87 256.73 299.5 326.72]

// A = matriz_mc(x,1)
// deter = det((A')*A)
// disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
[p1, erV, A, erN] = minimosCuadrados(x, y, 1)
disp("El polinomio de grado 1 es: ")
disp(p1)
disp("El error es: ")
disp(erN)
// A = matriz_mc(x,2)
// deter = det((A')*A)
// disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
[p2, erV, A, erN] = minimosCuadrados(x, y, 2)
disp("El polinomio de grado 2 es: ")
disp(p2)
disp("El error es: ")
disp(erN)
// A = matriz_mc(x,3)
// deter = det((A')*A)
// disp("La matriz A del método tiene rango completo, pues: det(A^T*A) = "+string(deter))
[p3, erV, A, erN] = minimosCuadrados(x, y, 3)
disp("El polinomio de grado 3 es: ")
disp(p3)
disp("El error es: ")
disp(erN)
rango = [min(x):0.0001:max(x)]
scf()
plot(rango, horner(p1, rango),"r")
plot(rango, horner(p2, rango),"g")
plot(rango, horner(p3, rango),"b")
plot(x', y', "r*")
h1 = legend(["Grado 1", "Grado 2", "Grado 3"])

// Ejercicio 9

/*
function w = f(x,n)
    nodos = -5:10/n:5
    y = 1./(1.+nodos.^2)
    pol = DD_Newton(nodos,y)
    aprox = horner(pol, x)
    w = 1./(1.+x.^2) - aprox
endfunction
*/
x = -5:0.01:5
//plot2d(x,[f(x,2)' f(x,4)' f(x,6)' f(x,10)' f(x,14)'],[2,3,4,5,6],leg="n=2@n=4@n=6@n=10@n=14")
for n=2:2:14
    if (n <> 8 && n <> 12)
        nodos = -5:10/n:5
        y = 1./(1.+x.^2)
        apy = 1./(1.+nodos.^2)
        pol = DD_Newton(nodos,apy)
        aprox = horner(pol, x)
        w = x - aprox
        scf()
        plot(x, aprox, "b")
        plot(x, y, "r")
        plot(x, w, "g")
        h1 = legend(["Aprox Grado "+ string(n), "1/1+x^2", "Error"])
    end
end

/*a=gca()
a.x_location = "origin"
a.y_location = "origin"
plot(x,f(x,2))
*/


/*
    Podemos ver que al aumentar el grado del polinomio, el error cerca 
    de 0 disminuye, pero en los extremos aumenta considerablemente.
    Considero que esto se debe al fenomeno de Runge.

*/

// Ejercicio 10

nodos = Cheb(4)'
rango = [-1:0.0001:1]
func = exp(rango)
pol = DD_Newton(nodos, exp(nodos))

scf()
plot(rango, horner(pol, rango), "r")
plot(rango, func, "g")
h1 = legend(['Polinomio Interpolacion','f(x)'])
a=get("current_axes")//get the handle of the newly created axes
a.axes_visible="on"; // makes the axes visible
a.y_location= "middle"
scf()
plot(rango, abs(func - horner(pol, rango)), "b")
h1 = legend(['Error del polinomio'])
a=get("current_axes")//get the handle of the newly created axes
a.axes_visible="on"; // makes the axes visible
a.y_location= "middle"


// Ejercicio 11

rango = [0:0.00001:%pi/2]
func = cos(rango)
grado = 3

nodos = NodosChebyshev(grado+1,0,%pi/2)

p = DD_Newton(nodos',cos(nodos'))

scf()
plot(rango,horner(p,rango),"r")
plot(rango,func,"b")
plot(rango,abs(horner(p,rango)-func),"g")
h1 = legend(["Aprox Grado 3", "cos(x)", "Error"])


// Casos Especiales

// Funcion exponencial natural. algunas veces conviene suponer que los datos tienen una relacion exponencial de la forma f(x) = a_2*e^(a_1*x)
// No es posible obtener una solucion exacta del sistema no lineal en a_1 y a_2.
// La solución es tomar el logaritmo de la funcion de aproximacion
// ln(f(x)) = ln(a_2) + a_1*x
// se puede definir el error como
// e_i =  ln(f(x_i)) - ln(y_i) = ln(a_2) + a_1*x_i -ln(y_i)
// e = Ax - b
// A = [1 x1; 1 x2; 1 x3; ... ; 1 xm]
// b = [ln(y1); ln(y2); ln(y3);...;ln(ym)]
// x = [ln(a_2); a_1]

// Funcion potencial
// f(x) = a_2*x^a_1
// ln(f(x)) = ln(a_2)+a_1 * ln(x)
// e_i = ln(f(x_i)) - ln(y_i) = ln(a_2) + a_1*ln(x_i) - ln(y_i)
// e = Ax - b
// A = [1 ln(x1); 1 ln(x2); 1 ln(x3); ... ; 1 ln(xm)]
// b = [ln(y1); ln(y2); ln(y3);...;ln(ym)]
// x = [ln(a_2); a_1]

// Factorizacion QR
// este esta mal? ayuda
function [Q, R] = factorizacionQR(A)
    [nA, mA] = size(A)
    v(1) = norm(A(:,1))
    Q(:,1) = A(:,1)/v(1)
    R = zeros(mA,mA)
    R(1,1) = v(1)
    for k = 2:mA
        suma = 0
        for i=1:k-1
            R(i,k) = A(:,k)'*Q(:,i)
            suma = suma + (A(:,k)'*Q(:, i))*Q(:, i)
        end
        comun = A(:,k) - suma
        v(k) = norm(comun)
        Q(:, k) = comun/v(k)
        R(k,k) = v(k)
    end
endfunction

function ayuda = minCuaQR(x, b, n)
    A = matriz_mc(x, n)
    [Q, R] = factorizacionQR(A)
    // superior(R, Q'*b')
    [ayuda,a] = gausselimPP(R,(Q')*(b'))
endfunction

function x = superior (A, b)
    [nA, mA] = size(A)
    [nb, mb] = size(b)
    
    if nA<>mA then
        error('superior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('superior - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    n = nA
    a = [A b]
    x(n) = a(n, n + 1) / a(n, n);
    for i = n-1:-1:1
        suma = 0
        for j = i+1:n
            suma = suma + a(i,j) * x(j)
        end
        x(i) = (a(i,n+1) - suma) / a(i, i)
    end
endfunction
