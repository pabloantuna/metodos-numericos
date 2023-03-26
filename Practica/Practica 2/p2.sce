function b = miHorner(p, x)
	n = degree(p);
	a = coeff(p);
	b = a(n);
	for i=1:(n-1)
		b = a(i) + x * b;
    end
endfunction

function [b, d] = miHornerEvolucionado(p,x)
	n = degree(p);
	a = coeff(p);
	b = a(n);
	if (n == 0) then
	d = 0;
	else
	d = a(n);
    end
	for i=1:(n-2)
		b = a(i) + x * b;
        d = b + x * d;
    end
    
    if (n > 1) then
        b = a(1) + x * b;
    end
endfunction

// funcion f es la ley de la función dada por un string, usa como 
// variable x
// v es el valor donde se evaluará la derivada
// n es el orden de derivación
// h es el paso de derivación
function d = derivar(f, v, n, h)
    deff("y=DF0(x)", "y="+f);
    if n == 0 then
        d = DF0(v);
    else
        for i=1:(n-1)
            funADerivar = "DF"+string(i-1);
            deff("y=DF"+string(i)+"(x)", "y=("+funADerivar+"(x+"+ string(h) + ")-" + funADerivar + "(x))/"+string(h));
        end
        ultimaFun = "DF"+string(n-1);
        deff("y=DFn(x)", "y=("+ultimaFun+"(x+"+ string(h) +")-" + ultimaFun + "(x))/" + string(h));
        d = DFn(v);
    end
endfunction

function d = derivarMasPiola(f, v, n, h)
    deff("y=DF0(x)", "y="+f);
    if n == 0 then
        d = DF0(v);
    else
        for i=1:(n-1)
            funADerivar = "DF"+string(i-1);
            deff("y=DF"+string(i)+"(x)", "y=numderivative("+funADerivar+",x,"+string(h)+"4)");
        end
        ultimaFun = "DF"+string(n-1);
        deff("y=DFn(x)", "y=numderivative("+ultimaFun+",x,"+ string(h) +"4)");
        d = DFn(v);
    end
endfunction

// Apartado a)

printf("Valor esperado para derivada de exponencial en 0: 1\n\n");
for i=1:4
    cociente = derivar("exp(x)",0,i,0.01);
    numder = derivarMasPiola("exp(x)",0,i,0.01)
    printf("Valor derivada orden %d con cociente: %.5e\n", i, cociente);
    printf("Valor derivada orden %d con numderivative: %.5e\n", i, numder);
    printf("Error relativo cociente orden %d: %.5e\n", i,(abs(cociente-1)/1));
    printf("Error relativo numderivate orden %d: %.5e\n", i,(abs(numder-1)/1));
    printf("\n");
end

esperado =  list(8+exp(1), 24+exp(1), 48+exp(1), 48+exp(1))
for i=1:4
    cociente = derivar("(2*x^4)+exp(x)",1,i,0.01)
    numder = derivarMasPiola("(2*x^4)+exp(x)",1,i,0.01)
    printf("Valor esperado para derivada orden %d de 2x^4 + exp(x) en 1: %.5e\n\n", i, esperado(i));
    printf("Valor derivada orden %d con cociente: %.5e\n", i, cociente);
    printf("Valor derivada orden %d con numderivative: %.5e\n", i, numder);
    printf("Error relativo cociente orden %d: %.5e\n", i, (abs(cociente-esperado(i))/esperado(i)));
    printf("Error relativo numderivate orden %d: %.5e\n", i, (abs(numder-esperado(i))/esperado(i)));
    printf("\n");
end

// Apartado b)

// En la implementación de cociente incremental, el error crece dependiendo del orden y del paso de derivación.
// Para un paso de derivación muy chico, en ordenes muy altos, se eleva a una potencia alta y por lo tanto
// cada vez se hace mas chico

function valor = taylor(f, v, n)
    h = 0.01

    for i = 1:n
        numerador = derivar(f, 0 , i, h)
        denominador = factorial(i)
        coeficientes(i) = numerador / denominador
    end
    p = poly(coeficientes, "x", "coeff")
    valor = horner(p, v) + f(0)
endfunction