//Ej 1
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

function x = inferior (A, b)
    [nA, mA] = size(A)
    [nb, mb] = size(b)
    
    if nA<>mA then
        error('inferior - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('inferior - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    n = nA
    a = [A b]
    x(1) = a(1, n + 1) / a(1, 1);
    for i = 2:n
        suma = 0
        for j = 1:(i-1)
            suma = suma + a(i,j) * x(j)
        end
        x(i) = (a(i,n+1) - suma) / a(i, i)
    end
endfunction


// Ejercicio 2

// Metodo dado por catedra GaussElim
function [x,a] = gausselim(A,b)
    // Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
    // dada la matriz de coeficientes A y el vector b.
    // La función implementa el método de Eliminación Gaussiana sin pivoteo.  
    
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A b]; // Matriz aumentada
    count = 0
    // Eliminación progresiva
    n = nA;
    for i=1:n-1
        for j=i+1:n
            mjk = a(j, i)/ a(i, i)
            a(j,i) = 0;
            a(j, (i+1):(n+mb)) = a(j, (i+1):(n+mb)) - mjk * a(i, (i+1):(n+mb))
            count = count + 1;
        end;
    end;
    
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    count = count + 1;
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
            count = count + 1
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
        count = count + 1;
    end;
    disp("La cantidad de operaciones fue: ");
    disp(count)

endfunction

// Ejemplos de aplicación
A = [3 -2 -1; 6 -2 2; -9 7 1]
b = [0 6 -1]'

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]'

[x,a] = gausselim(A,b)
disp(x)
disp(a)

A2 = [0 2 3; 2 0 3; 8 16 -1]
b2 = [7 13 -3]'

A2 = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b2 = [-8 -20 -2 4]'

//[x2,a2] = gausselim(A2,b2)

// !--error 27 
//Division by zero...
//at line      24 of function gausselim called by :  
//[x2,a2] = gausselim(A2,b2)

// apartado b
disp("Ejercicio 2b: ")
//b
//i
disp("Apartado i")
A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]'
[x,a] = gausselim(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

//ii
disp("Apartado ii")
A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]'
[x,a] = gausselim(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

//iii
disp("Apartado iii")
A = [1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2]
b = [2 1 0 -3]'
[x,a] = gausselim(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

//Ej 3
//Apartado a
// Eliminacion de Gauss sin Pivoteo
// Resuelve multiples sistemas de ecuaciones
function X = gauss_mult(A,B)
    [nA,mA] = size(A) 
    [nB,mB] = size(B)
    
    if nA<>mA then
        error('gaussmult - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nB then
        error('gaussmult - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A B]; // Matriz aumentada
    
    // Eliminación progresiva
    n = nA;
    for i = 1:(n-1)
        for j = (i+1):n
            mjk = a(j,i)/a(i,i)
            a(j,i)=0
            a(j,(i+1):(n+mB)) = a(j,(i+1):(n+mB)) - mjk*a(i,(i+1):(n+mB))
        end
    end
    X(n,1:mB) = a(n,(n+1):(n+mB))/a(n,n)
    for i = (nA-1):-1:1
      X(i,1:mB) = (a(i,(mA+1):(mA+mB)) - (a(i,(i+1):mA)*X((i+1):mA,1:mB)))/a(i,i)
    end 
endfunction

// apartado b
A = [1 2 3; 3 -2 1; 4 2 -1]
B = [14 9 -2; 2 -5 2; 5 19 12]

X = gauss_mult(A,B)
disp("La respuesta es")
disp(X)

// apartado c
// Calcula la inversa de una matriz.
function X = inversa(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('inversa - La matriz A debe ser cuadrada');
        abort;
    end
    I = eye(A)
    // la idea es aplicar el algoritmo
    // de eliminacion gaussiana para multiples sistemas
    // poniendo en lugar de un b, la matriz identidad
    X = gauss_mult(A,I)
endfunction

invA = inversa(A)
disp("La inversa es")
disp(invA)

// Ejercicio 4
function d = determinante(A)
    [n,m] = size(A) 
    if n<>m then
        error('determinante - La matriz A debe ser cuadrada');
        abort;
    end    
    // Eliminación progresiva
    d = 1
    for i = 1:n
        for j = (i+1):n
            mjk = A(j,i)/A(i,i)
            A(j,i)=0
            A(j,(i+1):n) = A(j,(i+1):n) - mjk*A(i,(i+1):n)
        end
        d = d * A(i,i)
    end
endfunction

// Ejercicio 5
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

// Ejemplo de aplicación
A2 = [0 2 3; 2 0 3; 8 16 -1]
b2 = [7 13 -3]'

[x2,a2] = gausselimPP(A2,b2)
disp(x2)
disp(a2)

// apartado b
disp("Ejercicio 5b: ")
//b
//i
disp("Apartado i")
A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1]
b = [4 1 -3 4]'
[x,a] = gausselimPP(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

//ii
disp("Apartado ii")
A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]'
[x,a] = gausselimPP(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

//iii
disp("Apartado iii")
A = [1 1 0 4; 2 1 -1 1; 4 -1 -2 2; 3 -1 -1 2]
b = [2 1 0 -3]'
[x,a] = gausselimPP(A,b)
disp("Solucion")
disp(x)
disp("Aumentada")
disp(a)

// Ejercicio 6

function x = tridiagonal(A, b)
    [nA, mA] = size(A)
    [nb, mb] = size(b)
    
    if nA<>mA then
        error('tridiagonal - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('tridiagonal - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A b] // aumentada
    n = nA
    for i=1:(n-1)
        mji = a(i+1,i)/a(i,i)
        a(i+1,i) = 0
        a(i+1,i+1) = a(i+1,i+1)-mji*a(i,i+1)
        a(i+1,n+1) = a(i+1,n+1)-mji*a(i,n+1)
    end
    //sust regres
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        x(i) = a(i,n+1)- a(i,i+1)*x(i+1)/a(i,i)
    end
endfunction

// Ejercicio 7
// Eliminacion de Gauss con Pivoteo Parcial para obtener la factorizacion PA = LU
function [L, U, P] = gaussLU(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('gaussLU - La matriz A debe ser cuadrada');
        abort;
    end
    
    I = eye(A)
    U = A
    L = I
    P = I
    n = nA
    for k=1:(n-1)
        kpivot = k; umax = abs(U(k,k));  //pivoteo
       for i=k+1:n
            if abs(U(i,k))>umax then
                kpivot = i; umax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:n); U(kpivot,k:n) = U(k,k:n); U(k,k:n) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        for j=k+1:n
            L(j,k) = U(j,k) / U(k,k)
            U(j, k:n) = U(j, k:n) - L(j,k) * U(k, k:n)
        end
    end
endfunction
A = [2 1 1 0; 4 3 3 1; 8 7 9 5; 6 7 9 8]
[L,U,P] = gaussLU(A)
disp("P")
disp(P)
disp("L")
disp(L)
disp("U")
disp(U)

// Ejercicio 8
// apartado a
A = [1.012 -2.132 3.104; -2.132 4.096 -7.013; 3.104 -7.013 0.014]
//mio
[L,U,P] = lu(A)
disp("APARTADO A")

disp("P mio")
disp(P)
disp("L mio")
disp(L)
disp("U mio")
disp(U)

//scilab
[L,U,P] = lu(A)

disp("P scilab")
disp(P)
disp("L scilab")
disp(L)
disp("U scilab")
disp(U)

// apartado b
A = [2.1756 4.0231 -2.1732 5.1967; -4.0231 6.0000 0 1.1973; -1.0000 5.2107 1.1111 0; 6.0235 7.0000 0 4.1561]
//mio
[L,U,P] = lu(A)
disp("APARTADO B")

disp("P mio")
disp(P)
disp("L mio")
disp(L)
disp("U mio")
disp(U)

//scilab
[L,U,P] = lu(A)

disp("P scilab")
disp(P)
disp("L scilab")
disp(L)
disp("U scilab")
disp(U)

// Ejercicio 9 en papel

// Ejercicio 10

function [L, U] = doolittle(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('doolittle - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(i, k) * U(k, j)
            end
            U(i,j) = A(i,j) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(m,k)*U(k,i)
                end
                L(m,i) = (A(m,i)-suma)/U(i,i)
            end
        end
    end
endfunction

// Resuelve un sistema de ecuaciones aplicando Doolittle
function [L,U,x]= d_solver(A, b)
    [L,U] = doolittle(A)
    y = inferior(L,b)
    x = superior(U,y)
endfunction

// aplicacion

A = [1 2 3 4; 1 4 9 16; 1 8 27 64; 1 16 81 256]
b = [2 10 44 190]'
[L, U, x] = d_solver(A, b)

disp("DOOLITTLE L")
disp(L)
disp("DOOLITTLE U")
disp(U)
disp("DOOLITTLE x")
disp(x)
// y = L\b
// x = U\y

// Ejercicio 11

function [U, ind] = CholeskyV1(A)
eps = 1.0e-8
n = size(A,1)
U = zeros(n,n)
for k = 1:n
    if k==1 then
        t = A(k,k)
    else 
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    end

    if t <= eps
        printf("Matriz no definida positiva.\n")
        ind = 0
        return
    end
    U(k,k)= sqrt(t)
    for j = k+1:n
        if k==1 then 
            U(k,j) = A(k,j)/U(k,k)
        else 
            U(k,j) = ( A(k,j) - U(1:k-1,k)' * U(1:k-1,j) )/U(k,k)
        end
    end
end
ind = 1
endfunction

A = [4 1 1; 1 2 2; 1 2 3]

[U,ind] = CholeskyV1(A)
disp("CHOLESKY V1")
disp(U)
disp(ind)

function [U,ind] = CholeskyV2(A)
// Factorización de Cholesky.
// Trabaja únicamente con la parte triangular superior.
//
// ind = 1  si se obtuvo la factorización de Cholesky.
//     = 0  si A no es definida positiva
//
//******************
eps = 1.0e-8
//******************

n = size(A,1)
U = zeros(n,n)

t = A(1,1)
if t <= eps then
    printf('Matriz no definida positiva.\n')
    ind = 0
    return
end
U(1,1) = sqrt(t)
for j = 2:n
    U(1,j) = A(1,j)/U(1,1)
end
    
for k = 2:n
    t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k)
    if t <= eps then
        printf('Matriz no definida positiva.\n')
        ind = 0
        return
    end
    U(k,k) = sqrt(t)
    for j = k+1:n
        U(k,j) = ( A(k,j) - U(1:k-1,k)'*U(1:k-1,j) )/U(k,k)
    end
end
ind = 1

endfunction

A = [4 1 1; 8 2 2; 1 2 3]
disp("CHOLESKY V2 a")
[U,ind] = CholeskyV2(A)
disp(U)
disp(ind)

disp("CHOLESKY V2 b")
B = [5 2 1 0; 2 5 2 0; 1 2 5 2; 0 0 2 5]
disp(B)
[U,ind] = CholeskyV2(B)
disp(U)
disp(ind)

disp("CHOLESKY V2 c")
C = [5 2 1 0; 2 -4 2 0; 1 2 2 2; 0 0 2 5]
disp(C)
[U,ind] = CholeskyV2(C)



// ejercicio 12
// apartado a

function x = chol_solverV1(A, b)
    [R, ind] = CholeskyV1(A)
    if ind==0 then
        error('chol_solverV1 - Matriz no definida positiva')
        abort;
    end
    Rt = R'
    g = inferior(Rt,b)
    x = superior(R,g)
endfunction

function x = chol_solverV2(A, b)
    [R, ind] = CholeskyV2(A)
    if ind==0 then
        error('chol_solverV2 - Matriz no definida positiva')
        abort;
    end
    Rt = R'
    g = inferior(Rt,b)
    x = superior(R,g)
endfunction

A = [16 -12 8; -12 18 -6; 8 -6 8]
b = [76 -66 46]'

disp("Cholesky solver V1")
x = chol_solverV1(A,b)
disp(x)
disp("Cholesky solver V2")
x = chol_solverV2(A,b)
disp(x)


// apartado b

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

function x = qr_solver(A,b)
    [Q,R] = factorizacionQR(A)
    x = superior(R, Q'*b)
endfunction

disp("FACTORIZACION QR")
//A = [1 -1 4 0; 1 4 -2 0; 1 4 2 0; 1 -1 0 0]
A = [16 -12 8; -12 18 -6; 8 -6 8]
b = [76 -66 46]'
//[Q,R] = factorizacionQR(A)
//disp("Q")
//disp(Q)
//disp("R")
//disp(R)
disp("Solucion")
x = qr_solver(A,b)
disp(x)

// Extras

// Resuelve un sistema de ecuaciones diagonal
// Recibe la matriz de coeficientes y el vector resultado 
function x = diagonal (A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('diagonal - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('diagonal - dimensiones incompatibles entre A y b');
        abort;
    end;
    x = b/diag(A)
endfunction

// Primero consigue la factorizacion LU con pivoteo
// Luego resuelve el sistema de ecuaciones triangulares
function [L, U, P, x] = gaussSolucionaLU(A,b)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('gaussSolucionaLU - La matriz A debe ser cuadrada');
        abort;
    end
    P = eye(A)
    L = eye(A)
    U = [A]
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(U(k,k));  //pivoteo
        for i=k+1:n
            if abs(U(i,k))>amax then
                kpivot = i; amax = U(i,k);
            end;
        end;
        temp = U(kpivot,k:mA); U(kpivot,k:mA) = U(k,k:mA); U(k,k:mA) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        
        for i=k+1:n
            L(i,k) = U(i,k)/U(k,k)
            U(i,k:mA) = U(i,k:mA)-L(i,k)*U(k,k:mA)
        end;
    end;
    c = P*b
    y = inferior(L, c)
    x = superior(U, y)
endfunction

function [L,U] = crout(A)
    [nA, mA] = size(A)
    if nA<>mA then
        error('crout - La matriz A debe ser cuadrada');
        abort;
    end
    
    L = eye(A)
    U = eye(A)
    for i=1:nA
        for j=i:nA
            suma = 0
            for k=1:i-1
                suma = suma + L(j,k)*U(k,i)
            end
            L(j,i)=A(j,i) - suma
            for m=i+1:nA
                suma = 0
                for k=1:i-1
                    suma = suma + L(i,k)*U(k,m)
                end
                U(i,m) = (A(i,m)-suma)/L(i,i)
            end
        end
    end
endfunction

function [L,U,x] = crout_solver(A,b)
    [L,U]=crout(A)
    y = inferior(L,b)
    x = superior(U,y)
endfunction
