// Metodo de la biseccion
// Se recibe una funcion continua en [a, b], los puntos a y b, y la tolerancia del error. no toma el maximo de iteraciones porque biseccion siempre converge

function r = biseccionIterativo(f, a, b, tol)
    c = (a+b)/2
    while b-c > tol

       if f(b)*f(c) <= 0 then
           a = c
       else
           b = c
       end 
      c = (a+b)/2
    end
    r = c


endfunction

function r = biseccionRecursivo(f,a, b, tol)
    c = (a + b) / 2
    if b - c <= tol then
        r = c
    else
        if f(b)*f(c) <= 0 then
           r = biseccion(f, c, b, tol, 0)
       else
           r = biseccion(f, a, c, tol, 0)
       end 
    end

endfunction

function y = f1(x)
    y = sin(x) - (x^2)/2
endfunction

function y = f2(x)
    y = exp(-x) - x^4
endfunction

function y = f3(x)
    y = log(x) - (x-1)
endfunction
xdata = linspace(-1,3,50)
ydata = f1(xdata)
plot2d(xdata,ydata)
e = 10^-2
printf("Metodo biseccion con item a: %f\n", biseccionIterativo(f1,-1,1,e))
printf("Metodo biseccion con item a: %f\n", biseccionIterativo(f1,1,2,e))
printf("Metodo biseccion con item b: %f\n", biseccionIterativo(f2,-1,1,e))
printf("Metodo biseccion con item c: %f\n", biseccionIterativo(f3,-1,1,e))
