// Metodo de Jacobi
function x = jacobi(A,b,x0,eps)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('jacobi - La matriz A debe ser cuadrada');
        abort;
    end
    for k=1:nA-1
        [v,i]=max(abs(A(k:nA,k)))
        kpivot = k-1+i
        temp = A(kpivot,:); A(kpivot,:) = A(k,:); A(k,:) = temp
        temp = b(kpivot,:); b(kpivot,:) = b(k,:); b(k,:) = temp
    end
    I = eye(nA,mA)
    N = diag(diag(A))
    T_j = I-inv(N)*A
    n(1) = norm(T_j, 1)
    n(2) = norm(T_j, 'inf')
    n(3) = norm(T_j, 'fro')
    n(4) = norm(T_j)
    if min(n) >= 1 then
        if max(abs(spec(T_j))) > = 1 then
            disp("La solucion no converge para todo punto inicial")
            x = %nan
            return
        end
    end
    
    n = nA
    y = x0
    x = zeros(n,1)
    iter = 0
    for k=1:n
        suma = 0
        for j=1:k-1
            suma = suma + A(k,j)*y(j)
        end
        for j=k+1:n
            suma = suma + A(k,j)*y(j)
        end
        x(k) =(b(k)-suma)/A(k,k)
    end
    iter = iter + 1
    while (norm(x-y) > eps)
        y = x
        for k=1:n
            suma = 0
            for j=1:k-1
            suma = suma + A(k,j)*y(j)
            end
            for j=k+1:n
                suma = suma + A(k,j)*y(j)
            end
            x(k) =(b(k)-suma)/A(k,k)
        end
        iter = iter + 1
    end
    printf("Iteraciones Jacobi: %d\n", iter)
endfunction

// Metodo de Gauss-Seidel
function x = gauss_seidel(A,b,x0,eps)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('jacobi - La matriz A debe ser cuadrada');
        return;
    end
    for k=1:nA-1
        [v,i]=max(abs(A(k:nA,k)))
        kpivot = k-1+i
        temp = A(kpivot,:); A(kpivot,:) = A(k,:); A(k,:) = temp
        temp = b(kpivot,:); b(kpivot,:) = b(k,:); b(k,:) = temp
    end
    I = eye(A)
    N = tril(A)
    T_gs = I-inv(N)*A
    n(1) = norm(T_gs, 1)
    n(2) = norm(T_gs, 'inf')
    n(3) = norm(T_gs, 'fro')
    n(4) = norm(T_gs)
    if min(n) >= 1 then
        if max(abs(spec(T_gs))) > = 1 then
            disp("La solucion no converge para todo punto inicial")
            x = %nan
            return
        end
    end
    
    
    n = nA
    x = x0
    y = x
    iter = 0
    for k=1:n
        suma = 0
        for j=1:k-1
            suma = suma + A(k,j)*x(j)
        end
        for j=k+1:n
            suma = suma + A(k,j)*x(j)
        end
        x(k) =(b(k)-suma)/A(k,k)
    end
    iter = iter + 1
    while (norm(x-y) > eps)
        y = x
        for k=1:n
            suma = 0
            for j=1:k-1
            suma = suma + A(k,j)*x(j)
            end
            for j=k+1:n
                suma = suma + A(k,j)*x(j)
            end
            x(k) =(b(k)-suma)/A(k,k)
        end
        iter = iter + 1
    end
    printf("Iteraciones Gauss-Seidel: %d\n", iter)
endfunction

// Chequea si una matriz es diagonal dominante
// devuelve 1 si es asi 0 en caso contrario
function x = diagonalDominante(A)
    [nA,mA] = size(A)
    for i = 1:nA
        suma = 0
        for j = 1:i-1
            suma = suma + abs(A(i,j))
        end
        for j = i+1:mA
            suma = suma + abs(A(i,j))
        end
        if suma >= abs(A(i,i))
            x = 0
            return
        end
    end
    x = 1
endfunction

// Ej 3

function M = matrizGaussSeidel(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('jacobi - La matriz A debe ser cuadrada');
        return;
    end

    I = eye(A)
    N = tril(A)
    M = I-inv(N)*A
endfunction

// a es el valor para la diag inf
// b para la diag media
// c para la diag sup
// n el tamaño de la matriz
function M = construirTriDiagonal(a, b, c, n)
    M = zeros(n, n)
    M(1,1)=b
    M(1,2)=c
    for i = 2:n-1
        M(i,i-1) = a
        M(i,i) = b
        M(i,i+1) = c
    end
    M(n,n-1) = a
    M(n,n) = b
endfunction

A = construirTriDiagonal(-1, 2, -1, 5)
M = matrizGaussSeidel(A)
disp("ACA")
disp(A)
disp(M)
A = construirTriDiagonal(-1, 2, -1, 8)
M = matrizGaussSeidel(A)
disp("ACA")
disp(A)
disp(M)
A = construirTriDiagonal(-1, 2, -1, 10)
M = matrizGaussSeidel(A)
disp("ACA")
disp(A)
disp(M)
A = construirTriDiagonal(-1, 2, -1, 15)
M = matrizGaussSeidel(A)
disp("ACA")
disp(A)
disp(M)
A = construirTriDiagonal(-1, 2, -1, 20)
M = matrizGaussSeidel(A)
disp("ACA")
disp(A)
disp(M)

//ej 4
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

N = 100
A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3)
b = ones(N,1)

disp("N = 100")

disp("GAUSS ELIM PP")
tic(); gausselimPP(A,b); t = toc()
disp(t)

disp("GAUSS SEIDEL eps 10^-6")
tic(); gauss_seidel(A,b,zeros(N,1),0.000001); t = toc()
disp(t)

disp("GAUSS SEIDEL eps 10^-11")
tic(); gauss_seidel(A,b,zeros(N,1),0.00000000001); t = toc()
disp(t)

N = 500
A = 8*eye(N,N) + 2*diag(ones(N-1,1),1) + 2*diag(ones(N-1,1),-1) + diag(ones(N-3,1),3) + diag(ones(N-3,1),-3)
b = ones(N,1)

disp("N = 500")

disp("GAUSS ELIM PP")
tic(); gausselimPP(A,b); t = toc()
disp(t)

disp("GAUSS SEIDEL eps 10^-6")
tic(); gauss_seidel(A,b,zeros(N,1),0.000001); t = toc()
disp(t)

disp("GAUSS SEIDEL eps 10^-11")
tic(); gauss_seidel(A,b,zeros(N,1),0.00000000001); t = toc()
disp(t)

// Metodo de Gauss-Seidel con relajacion
function x = gauss_relajacion(A,b,x0,w,eps)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('gauss-seidel relajacion - La matriz A debe ser cuadrada');
        abort;
    end
    for k=1:nA-1
        [v,i]=max(abs(A(k:nA,k)))
        kpivot = k-1+i
        temp = A(kpivot,:); A(kpivot,:) = A(k,:); A(k,:) = temp
        temp = b(kpivot,:); b(kpivot,:) = b(k,:); b(k,:) = temp
    end
    I = eye(nA,mA)
    N = A
    for i = 1:nA-1
        for j = i+1:nA
            N(i,j) = 0
        end
    end
    inversaN = inv(N)
    Norma = I-inversaN*A
    n(1) = norm(Norma, 1)
    n(2) = norm(Norma, 'inf')
    n(3) = norm(Norma, 'fro')
    n(4) = norm(Norma)
    if min(n) >= 1 then
        if max(abs(spec(Norma))) > = 1 then
            disp("La solucion no converge para todo punto inicial")
            x = %nan
            abort
        end
    end
    
    n = nA
    x = x0
    y = x
    iter = 0
    for k=1:n
        suma = 0
        for j=1:k-1
            suma = suma + A(k,j)*x(j)
        end
        for j=k+1:n
            suma = suma + A(k,j)*x(j)
        end
        x(k) =(1-w)*x(k)+ w*(b(k)-suma)/A(k,k)
    end
    iter = iter + 1
    while (norm(x-y) > eps)
        y = x
        for k=1:n
            suma = 0
            for j=1:k-1
            suma = suma + A(k,j)*x(j)
            end
            for j=k+1:n
                suma = suma + A(k,j)*x(j)
            end
            x(k) =(1-w)*x(k)+ w*(b(k)-suma)/A(k,k)
        end
        iter = iter + 1
    end
    printf("Iteraciones Gauss-Seidel Relajacion: %d\n", iter)
endfunction

// basicamente lo mismo que calcular el w con la formula
// para tridiag pero reutilizando matrices para hacer calculos una sola vez
function x = gauss_relajacion_tri(A,b,x0,eps)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('gauss-seidel relajacion - La matriz A debe ser cuadrada');
        abort;
    end
    
    I = eye(nA,mA)
    N = A
    for i = 1:nA-1
        for j = i+1:nA
            N(i,j) = 0
        end
    end
    inversaN = inv(N)
    Norma = I-inversaN*A
    n(1) = norm(Norma, 1)
    n(2) = norm(Norma, 'inf')
    n(3) = norm(Norma, 'fro')
    n(4) = norm(Norma)
    if min(n) >= 1 then
        if max(abs(spec(Norma))) > = 1 then
            disp("La solucion no converge para todo punto inicial")
            x = %nan
            abort
        end
    end
    D = diag(diag(A))
    T = I-inv(D)*A
    p = max(abs(spec(T)))
    w = 2/(1+sqrt(1-p^2))
    
    n = nA
    x = x0
    y = x
    iter = 0

    x(1) =(1-w)*x(1)+ w*(b(1)-A(1,2)*x(2))/A(1,1)
    for k=2:n-1
        suma = A(k,k-1)*x(k-1)+ A(k,k+1)*x(k+1)
        x(k) =(1-w)*x(k)+ w*(b(k)-suma)/A(k,k)
    end
    x(n) =(1-w)*x(n)+ w*(b(n)-A(n,n-1)*x(n-1))/A(n,n)
    iter = iter + 1
    while (norm(x-y) > eps)
        y = x
        x(1) =(1-w)*x(1)+ w*(b(1)-A(1,2)*x(2))/A(1,1)
        for k=2:n-1
            suma = A(k,k-1)*x(k-1)+ A(k,k+1)*x(k+1)
            x(k) =(1-w)*x(k)+ w*(b(k)-suma)/A(k,k)
        end
        x(n) =(1-w)*x(n)+ w*(b(n)-A(n,n-1)*x(n-1))/A(n,n)
        iter = iter + 1
    end
    printf("Iteraciones Gauss-Seidel Relajacion Tridiagonal: %d\n", iter)
endfunction

