function [radios, centros] = cotas(A)
    [n, m] = size(A)
    centros = diag(A)
    radios = sum(abs(A), 'c') - abs(centros) // sumo todos los valores de la fila y le saco el valor abs de la diag (para no contarlo a el mismo)
    
    for i=1:n
        printf("El autovalor %d puede encontrarse en el circulo de Gerschgorin de radio %f y centro %f\n", i, radios(i), centros(i))
    end
endfunction

function circ(r,x,y)
     rect = [x-r,y-r,x+r,y+r]
    p = 0
    xgrid(1)
    plot2d(p,rect=rect)
    xarc(x-r,y+r,r*2,r*2,0,360*64)
    c = 0
//    xarc(x-r,y+r,2*r,2*r,0,360*64)
endfunction

function circa(r,x,y)
    xarc(x-r,y+r,2*r,2*r,0,360*64)
endfunction

function Gers(A)
    [n,m] = size(A);
    centros = diag(A);
    radios = sum(abs(A),'c') - abs(centros) ;
    
    // buscamos calcular un rectángulo que contenga a todos los circulos
    // esquina inferiro izquierda
    
    mx = round (min(centros - radios)-1);
    my = round (min(-radios)-1);
    
    // esquina superior derecha
    
    Mx = round(max(centros+radios)+1);
    My = round(max(radios)+1);
    
    rectangulo = [mx my Mx My];
    
    // dibujamos los autovalores
    plot2d(0,0,-1,"031","",rectangulo)
    replot(rectangulo); // reeplaza al rect
    xgrid();
    
    for i=1:n
        circ(radios(i),centros(i),0)
    end
    
endfunction

function CircGersValor(A)
    [n,m] = size(A);
    centros = diag(A);
    radios = sum(abs(A),'c') - abs(centros) ;
    
    // buscamos calcular un rectángulo que contenga a todos los circulos
    // esquina inferiro izquierda
    
    mx = round (min(centros - radios)-1);
    my = round (min(-radios)-1);
    
    // esquina superior derecha
    
    Mx = round(max(centros+radios)+1);
    My = round(max(radios)+1);
    
    rectangulo = [mx my Mx My];
    
    // dibujamos los autovalores
    plot2d(real(spec(A)),imag(spec(A)),-3,"031","",rectangulo)
    replot(rectangulo); // reeplaza al rect
    xgrid(4, 1.5, 9);
    
    for i=1:n
        circa(radios(i),centros(i),0)
    end
    
endfunction


// ejercicio 1
function ej1()
    A = [1 0 0; -1 0 1; -1 -1 2]
    B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
    C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
    D = [4 -1 0; -1 4 1; -1 -1 4]
    E = [3 2 1; 2 3 0; 1 0 3]
    F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]
    
    // A
    printf ("Matriz A\n\n")
    cotas(A)
    eigenS = spec(A)
    printf("Autovalores de scilab:")
    disp(eigenS)
    
    // B
    printf ("Matriz B\n\n")
    cotas(B)
    eigenS = spec(B)
    printf("Autovalores de scilab:")
    disp(eigenS)
    
    // C
    printf ("Matriz C\n\n")
    cotas(C)
    eigenS = spec(C)
    printf("Autovalores de scilab:")
    disp(eigenS)
    
    // D
    printf ("Matriz D\n\n")
    cotas(D)
    eigenS = spec(D)
    printf("Autovalores de scilab:")
    disp(eigenS)
    
    // E
    printf ("Matriz E\n\n")
    cotas(E)
    eigenS = spec(E)
    printf("Autovalores de scilab:")
    disp(eigenS)
    
    // F
    printf ("Matriz F\n\n")
    cotas(F)
    eigenS = spec(F)
    printf("Autovalores de scilab:")
    disp(eigenS)
endfunction



// este es el dado por la catedra
function gres(A)
    [n,m] = size(A);
    centros = diag(A);
    radios = sum(abs(A),'c') - abs(centros) ;
    
    // buscamos calcular un rectángulo que contenga a todos los circulos
    // esquina inferiro izquierda
    
    mx = round (min(centros - radios)-1);
    my = round (min(-radios)-1);
    
    // esquina superior derecha
    
    Mx = round(max(centros+radios)+1);
    My = round(max(radios)+1);
    
    rectangulo = [mx my Mx My];
    
    // dibujamos los autovalores
    plot2d(real(spec(A)),imag(spec(A)),-1,"031","",rectangulo)
    replot(rectangulo); // reeplaza al rect
    xgrid();
    
    for i=1:n
        circ(radios(i),centros(i),0)
    end
    
endfunction

//ejercicio 3

function A = crearMatrizAeps(eps)
    A = [1 -1 0; -2 4 -2; 0 -1 1+eps]
endfunction

function ej3()
    eps = 0.1
    for k=0:10
        epsR = eps*k
        printf("Epsilon = %f\n", epsR)
        A = crearMatrizAeps(epsR)
        p = poly(A,"x")
        printf("Polinomio caracteristico:")
        disp(p)
        printf("Raices del polinomio caracteristico:")
        disp(roots(p))
        eigenS = spec(A)
        printf("Autovalores de A:")
        disp(eigenS)
    end
endfunction

//ejercicio 5

function [x, z, iters] = metodoPotenciaRaro(A, z0, eps, maxiter)
    //iteracion 1
    w1 = A*z0
    z1 = w1 / norm(w1, 'inf')
    [wk, k] = max(abs(w1))      // si el maximo de absoluto me da 0 entonces w1 es el vector nulo.
    if z0(k) == 0 then
        k = 1
        while z0(k) == 0 && k <= nw
            k = k+1
        end
    end

    wk = w1(k)
    l1 = wk / z0(k)
    
    //iteracion 2
    w2 = A*z1
    z2 = w2 / norm(w2, 'inf')
    [wk, k] = max(abs(w2))
    if z0(k) == 0 then
        k = 1
        while z0(k) == 0 && k <= nw
            k = k+1
        end
    end
    wk = w2(k)
    l2 = wk / z1(k)
    
    //iteracion 3
    w3 = A*z2
    z3 = w3 / norm(w3, 'inf')
    [wk, k] = max(abs(w3))
    if z0(k) == 0 then
        k = 1
        while z0(k) == 0 && k <= nw
            k = k+1
        end
    end
    wk = w3(k)
    l3 = wk / z2(k)
    
    i = 3
    errN = (l3 - l2) / (l2-l1)
    errD = (1-errN)*(l3-l2)
    err = errN / errD
    // resto iteraciones
    while (i <= maxiter && err > eps)
        z2 = z3
        w3 = A*z3
        z3 = w3 / norm(w3, 'inf')
        l1 = l2
        l2 = l3
        [wk, k] = max(abs(w3))
        if z0(k) == 0 then
            k = 1
            while z0(k) == 0 && k <= nw
                k = k+1
            end
        end
        wk = w3(k)
        l3 = wk / z2(k)
        errN = (l3 - l2) / (l2-l1)
        errD = (1-errN)*(l3-l2)
        err = errN / errD
        i = i+1
    end
    x = l3
    z = z3
    iters = i
endfunction

function [v, zn, iters]= metodoPotencia(A, z0, eps, maxiter)
    V = max(abs(spec(A)))
    v = 0
    iters = 1
    w = A*z0
    zn = w/norm(w)
    [m,j] = max(abs(w))
    [nw, mw] = size(w)
    if z0(j) == 0 then
        j = 1
        while z0(j) == 0 && j <= nw
            j = j+1
        end
    end
    m = w(j)
    v = m/z0(j) 
    er1 = abs(V-v)
    er2 = norm(zn-z0)
    err = max(er1,er2)
    while (iters <= maxiter && err > eps)
        z0 = zn
        w = A*z0
        zn = w/norm(w,'inf')
        [m,j] = max(abs(w))
        if z0(j) == 0 then
            j = 1
            while z0(j) == 0 && j <= nw
                j = j+1
            end
        end
        m = w(j)
        v = m/z0(j)
        er1 = abs(V-v)
        er2 = norm(zn-z0,'inf')
        err = max(er1,er2)
        iters = iters+1
    end
endfunction

// aplicado
function ej5v1()
    A1 = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
    A2 = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
    Z0 = ones(4,1)
    eps = 10^(-8)
    maxiter = 1000
    [x, z, iters] = metodoPotenciaRaro(A1, Z0, eps, maxiter)
    printf("Matriz A1\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    printf("\n")
    
    [x,z, iters] = metodoPotenciaRaro(A2,Z0,eps,maxiter)
    printf("Matriz A2\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
endfunction

function ej5v2()
    A1 = [6 4 4 1; 4 6 1 4; 4 1 6 4; 1 4 4 6]
    A2 = [12 1 3 4; 1 -3 1 5; 3 1 6 -2; 4 5 -2 -1]
    Z0 = ones(4,1)
    eps = 10^(-8)
    maxiter = 1000
    [x, z, iters] = metodoPotencia(A1, Z0, eps, maxiter)
    printf("Matriz A1\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    printf("\n")
    
    [x,z, iters] = metodoPotencia(A2,Z0,eps,maxiter)
    printf("Matriz A2\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
endfunction

function papopepo()
    A = [1 0 0; -1 0 1; -1 -1 2]
    B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
    C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
    D = [4 -1 0; -1 4 1; -1 -1 4]
    E = [3 2 1; 2 3 0; 1 0 3]
    F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]
    
    Z0 = ones(3,1)
    eps = 10^(-8)
    maxiter = 1000
    [x, z, iters] = metodoPotencia(A, Z0, eps, maxiter)
    printf("Matriz A1\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    printf("\n")
    
    [x,z, iters] = metodoPotencia(B,Z0,eps,maxiter)
    printf("Matriz B\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotencia(C,Z0,eps,maxiter)
    printf("Matriz C\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotencia(D,Z0,eps,maxiter)
    printf("Matriz D\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotencia(E,Z0,eps,maxiter)
    printf("Matriz E\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotencia(F,Z0,eps,maxiter)
    printf("Matriz F\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
endfunction

function papopepo2()
    A = [1 0 0; -1 0 1; -1 -1 2]
    B = [1 0 0; -0.1 0 0.1; -0.1 -0.1 2]
    C = [1 0 0; -0.25 0 0.25; -0.25 -0.25 2]
    D = [4 -1 0; -1 4 1; -1 -1 4]
    E = [3 2 1; 2 3 0; 1 0 3]
    F = [4.75 2.25 -0.25; 2.25 4.75 1.25; -0.25 1.25 4.75]
    
    Z0 = ones(3,1)
    eps = 10^(-8)
    maxiter = 1000
    [x, z, iters] = metodoPotenciaRaro(A, Z0, eps, maxiter)
    printf("Matriz A1\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    printf("\n")
    
    [x,z, iters] = metodoPotenciaRaro(B,Z0,eps,maxiter)
    printf("Matriz B\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotenciaRaro(C,Z0,eps,maxiter)
    printf("Matriz C\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotenciaRaro(D,Z0,eps,maxiter)
    printf("Matriz D\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotenciaRaro(E,Z0,eps,maxiter)
    printf("Matriz E\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
    
    [x,z, iters] = metodoPotenciaRaro(F,Z0,eps,maxiter)
    printf("Matriz F\n")
    printf("Iteraciones realizadas: %d\n", iters)
    printf("Autovalor dominante: %f\n", x)
    printf("Autovector asociado:")
    disp(z)
endfunction
