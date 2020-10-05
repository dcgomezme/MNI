### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ 20e7ed8e-0593-11eb-1afb-7d77519ef9ee
begin
	using Images
	viga_path = "img/viga.png";viga = load(viga_path)
end

# ╔═╡ bbe25158-05ac-11eb-3cb6-ed5e87450d02
begin
	using Plots
	y = 0:0.1:2;
	plot(y,(y -> begin ((3*y+((y^2)/2))^3/(3+y)) - 40.7747 end),framestyle=:zerolines)
end

# ╔═╡ 7dbdf4d8-0587-11eb-0ce5-d397912d3eee
md"""
# Raíces y sistemas de Ecuaciones
"""

# ╔═╡ 119fe18e-0588-11eb-088b-8120e27e9433
md"""
Para cada uno de los ejercicios que se presentan a continuación, usted debe explicar el fenómeno y determinar las ecuaciones que gobiernan el problema. Luego, utilizar programas en computador para resolver lo que se pide. La entrega debe ser un único documento en formato pdf con la solución de cada ejercicio más un script de Octave/Matlab para cada ejercicio (4 en total). Tenga en cuenta que puede escribir funciones de ayuda para encontrar la solución a los ejercicios. En dicho caso,debe incluir estos archivos en su entrega.
"""

# ╔═╡ 648d66d2-0588-11eb-2b95-4517fc61e165
md"""
## I Viga
Con base en las dimensiones y cargas que es especifican en la siguiente figura, use el método de la bisección para encontrar la posición dentro de la viga donde no hay momento. Asuma la tolerancia que considere apropiada e indique el número de iteraciones que necesitó para resolver el ejercicio.
"""

# ╔═╡ 569072c0-0594-11eb-05e2-1985663cc0cc
md"""
### Desarrollo
Resolviendo la viga encontramos que las ecuaciones de momento flector que gobiernan el comportamiento en la viga son las siguientes:

|Tramo| M(x) |
|:----|-------:|
|$0\leq x < 3$|$\displaystyle 265x-\frac{50x^3}{9}$|
|$3\leq x < 6$|$\displaystyle -50x^2+415x-150$|
|$6\leq x <10$|$\displaystyle -185x + 1650$|
|$10\leq x\leq12$|$\displaystyle 100x - 1200$|

Si graficamos el diagrama de Momento flector obtendremos:
"""

# ╔═╡ 0c4f5008-0596-11eb-1db3-c14571aee5bb
flector_path = "img/flector.png";flector = load(flector_path)	

# ╔═╡ 381b8ac4-0598-11eb-23bc-03663f91c2f3
md""" 
Por inspección vemos que el lugar donde el momento flector es igual a cero se encuentra en el tramos $ C \leq x < D $ el cuál equivale al tramo $ 6 \leq x < 10$ , por lo tanto la ecuación que vamos a usar para determinar el valor de $x$ donde el momento flector es igual a $0$ será:

$$M(x) = -185x+1650$$
"""

# ╔═╡ aa3d2694-0598-11eb-3b34-11a65f96e210
function  biseccion(xl,xu,f,tol,ea=1)
	iter = 0;xra = xu;
	while ea > tol
		xr = (xl+xu)/2;
		ea = abs((xr - xra)/xr);
		if f(xl)*f(xr) < 0
			xu = xr;
		end
		if f(xr)*f(xu)< 0
			xl = xr;
		end
		xra = xr;
		iter += 1;
	end
	return xra, iter
end

# ╔═╡ 76a9ce70-059a-11eb-13b5-13c94ec1c3e0
begin
	M(x) = -185*x +1650;
	xra, iter = biseccion(6,10,M,0.01);
	"A una distancia de $xra pies el momento flector es igual a 0, el método le tomo $iter iteraciones encontrar este valor"
end

# ╔═╡ cc74b18a-059a-11eb-3167-8525aeee8f3f
md"""
## II Canal
Por un canal trapezoidal fluye agua a una tasa de $Q = 20m^3/s$. La profundidad crítica $ y $ para dicho canal satiface la ecuación.

$$0 = 1 - \frac{Q^2}{gA_c^3}B$$

Resuelva para la profundidad crítica con el uso de los métodos.

	a. Gráfico
	b. Bisección (con xl = 0.5 y xu = 2.5)
	c. Newton-Raphson (con x0 =0.5)

Haga el ejercicio hasta que el error aproximado caiga por debajo del 1 % o el número de iteraciones supere a 10. Analice sus resultados
"""

# ╔═╡ 684870e0-059d-11eb-3df9-21dba1e9dce7
md"""

### Desarrollo
Lo primero que vamos a hacer es encontrar una expresión más manejable para poder realizar los cálculos de una mejor manera.

$$0 = 1 -\frac{Q^2}{gA_c^3}B  \hspace{3mm}\rightarrow\hspace{3mm} \frac{Q^2}{gA_c^3}B = 1 \hspace{3mm}\rightarrow\hspace{3mm} \frac{Q^2}{g} = \frac{A_c^3}{B} \hspace{3mm}\rightarrow\hspace{3mm} 0 = \frac{A_c^3}{B} - \frac{Q^2}{g}$$

Si realizamos el reemplazo de las expresiones $\displaystyle A_c = 3y + \frac{y^2}{2} \hspace{3mm} \& \hspace{3mm} B = 3 + y$  la función a evaluar sería:

$$f_y = 0 = \frac{\big(3y+\frac{y^2}{2}\big)^3}{3 + y} - \frac{20^2}{9.81} \hspace{3mm}\rightarrow\hspace{3mm} f_y = \frac{\big(3y+\frac{y^2}{2}\big)^3}{3 + y} - 40.7747$$

Con el fin de aplicar el método de Newton - Raphson calculamos la derivada de la función:

$$f'_y = \frac{(3y+9)\big(\frac{y^2}{2}+3y\big)^2}{y+3} - \frac{\big(\frac{y^2}{2} + 3y \big)^3}{(y+3)^2}$$

#### a) Graficamente

Primero graficamos la función para poder encontrar el punto aproximado y graficamente ver donde se encuentra la raiz

Donde podemos ver que la raiz se encuentra aproximadamente en el valor de $y = 1.51$
"""

# ╔═╡ 583ec42c-05ae-11eb-1a75-051d43500d23
md"""
### b) Bisección
Teniendo formulada la función solamente debemos cambiar la función a evaluar y los límites del intervalo.
"""

# ╔═╡ 095a0834-05af-11eb-2a36-d1d1e179e27d
begin
	Cb(y) = ((3*y+((y^2)/2))^3/(3+y)) - 40.7747;
	yra , numiter = biseccion(0.5,2.5,Cb,0.01);
	"La profundidad crítica encontrada con el método de la bisección fue $yra luego de $numiter iteraciones"
end

# ╔═╡ f25fb776-05b0-11eb-37bd-09b89dc26a2f
md"""
### c) Newton - Raphson
Primero se debe formular el método dentro de una función y luego llamar la función para implementar el método en nuestro caso específico
"""

# ╔═╡ 5291f79e-05b1-11eb-0233-05ca4a69411b
function Newton_Raphson(x0,f,df,tol,imax)
	var =0;iter = 0;
	for i = 1 : imax;
		xi = x0 - (f(x0)/df(x0));
		ea = abs((xi - x0)/xi);
		if ea < tol; var = 1; break end
		x0 = xi; iter += 1;
	end
	if var == 0; "La función no converge" end
	return x0, iter
end

# ╔═╡ 4c4d09fe-05b2-11eb-1609-b1e4adb4b50e
begin
	Cn(y) = ((3*y+((y^2)/2))^3/(3+y)) - 40.7747;
	DCn(y) = (3*y + 9)*(y^2/2 + 3*y)^2/(y + 3) - (y^2/2 + 3*y)^3/(y + 3)^2;
	yc, itera = Newton_Raphson(0.5,Cn,DCn,0.01,10);
	"La profundidad crítica encontrada con el método de Newton-Raphson fue $yc luego de  $itera iteraciones"
end

# ╔═╡ a9445d72-05b4-11eb-0c50-e7b437caf6f2
md"""
## III Tubos
Un fluido se bombea en la red de tubos que se muestra en la figura a continuación.
"""

# ╔═╡ 0ab3162a-05b5-11eb-0cc6-1bce617dfa96
tubos_path = "img/tubos.png";tubos = load(tubos_path)

# ╔═╡ 29365c30-05b6-11eb-331c-19a485128ce2
md"""
En estado estacionario, se cumplen los balances de flujo siguientes:

$$Q_1 = Q_2 + Q_3$$
$$Q_3 = Q_4 + Q_5$$
$$Q_5 = Q_6 + Q_7$$

donde $Q_i$ es el flujo en el tubo $i[m^3/s]$. Además, la caída de presión alrededor de los tres lazos en los que el flujo es hacia la derecha debe ser ser igual a cero. La caída de presión en cada tramo de tubo circular se calcula por medio de la ecuación:

$$\Delta P = \frac{16}{\pi^2}\frac{fL\rho}{2D^5}Q^2$$

donde $\Delta P = $ caída de presión [Pa], $f =$ factor de fricción [-], $L =$ longitud del tubo [m], $\rho =$ densidad del fluido [kg/m$^3$], y $D =$ diámetro del tubo [m]. Escriba un programa (o desarrolle un algortmo en un *script*) que permita calcular el flujo en cada tramo del tubo, dado que $Q_1 = 1m^3/s$ y $\rho = 1.23kg/m^3$. Todos los tubos tienen $D = 500mm$ y $f = 0.005$. Las longitudes de los tubos son: $L_3 = L_5 = L_8 = L_9 = 2m; L_2 = L_4 = L_6 = 4m;$ y $L_7 = 8m$.

"""

# ╔═╡ 6c699394-065a-11eb-26d7-1b895b06b2d4
md"""
### Desarrollo

Gracias a los datos que nos da el problema es fácil ver que las incognitas a encontrar son la cantidad de caudales dentro del sistema menos 1, pues el valor de $Q_1 = 1m^3/s$, esto nos da un numero de incognitas igual a 9, por lo tanto debemos tener un sistema de 9 ecuaciones lineales con nueves incognitas. El método lineal nos permite convertir la ecuacion no lineal dada para $\Delta P$ en una ecuación lineal, asi que usaremos este método para resolver el sistema.
"""

# ╔═╡ 27dd8786-065e-11eb-0f4b-1ddcb3e6c6f9
md"""
#### Ecuaciones de Circuito

$$F(𝐈) = -C_2(Q_2^2) + C_3(Q_3^2) + C_4(Q_4^2) + C_9(Q_9^2) = 0$$
$$F(𝐈𝐈) = -C_4(Q_4^2) + C_5(Q_5^2) + C_6(Q_6^2) + C_8(Q_9^2) = 0$$
$$F(𝐈𝐈𝐈) = -C_6(Q_6^2) + C_7(Q_7^2) = 0$$

Donde:

$$\displaystyle C_i = \frac{16}{\pi^2}\frac{fL_i\rho}{2D^5}$$

Para linealizar las ecuaciones cada término $\pm C_i(Q_i^2)$ lo convertiremos en el término $\pm C_iQ_i^*Q_i$, donde el valor de $Q_i^*$ se le conoce como el canal supuesto el cual tiene un valor semilla arbitrario, pero n-ésima iteración su valor será igual a:

$$Q_n = \frac{Q_{n-1} - Q_{n-2}}{2}$$ 

Es decir el promedio entre el Caudal asumido en la iteración anterior y el Caudal obtenido en la iteración anterior.
"""

# ╔═╡ 0fef8ab6-065e-11eb-3903-cfbc89172dad
md"""
### Ecuaciones de nodo

Partimos de tres ecuaciones que el problema nos da y luego planteamos otras 3 para sumar 9 ecuaciones en total, 6 de nodo y 3 de circuito.

$$Q_1 = Q_2 + Q_3$$
$$Q_3 = Q_4 + Q_5 \hspace{3mm}\rightarrow\hspace{3mm} 0 = -Q_3 + Q_4 + Q_5$$
$$Q_5 = Q_6 + Q_7 \hspace{3mm}\rightarrow\hspace{3mm} 0 = -Q_5 + Q_6 + Q_7$$
$$Q_8 = Q_6 + Q_7 \hspace{3mm}\rightarrow\hspace{3mm} 0 =  Q_6 + Q_7 - Q_8$$
$$Q_9 = Q_4 + Q_8 \hspace{3mm}\rightarrow\hspace{3mm} 0 =  Q_4 + Q_8 - Q_9$$
$$Q_{10}= Q_2 + Q_9 \hspace{3mm}\rightarrow\hspace{3mm} 0 =  Q_2 + Q_9 - Q_{10}$$

"""



# ╔═╡ 876e567e-0660-11eb-2906-91173842bfad
md"""
### Sistema Matricial

$$\begin{bmatrix}
-C_2Q_2^* & C_3Q_3^* & C_4Q_4^* & 0 & 0 & 0 & 0 & C_9Q_9^*&0  \\
0 & 0 & -C_4Q_4^* & C_5Q_5^* & C_6Q_6^* & 0 & C_8Q_8^* & 0 & 0 \\
0 & 0 & 0 & 0 & -C_6Q_6^* & C_7Q_7^* & 0 & 0 & 0 \\
1 & 1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & -1 & 1 & 1 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & -1 & 1 & 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 1 & -1 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 & 1 & -1 & 0 \\
1 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & -1 \\
\end{bmatrix}
\begin{bmatrix} Q_2\\Q_3\\Q_4\\Q_5\\Q_6\\Q_7\\Q_8\\Q_9\\Q_{10}\\ \end{bmatrix}=
\begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ \end{bmatrix}$$

El valor de 1 en el vector de términos independientes es iguál al caudal dado de $Q_1 = 1m^3/s$ por lo tanto estas serán las unidades de los caudales obtenidos.

A continuación se planterá la función que resolverá el sistema matricial, y luego se planteará el sistema y se solucionará con la función escrita.
"""

# ╔═╡ aa576bf0-0663-11eb-1b89-35d73b118f51
function Sol(M,b,Q0,tol,r,imax,ea=1)
Q = Q0*ones(1,size(M)[2]);
  for i = 1 : imax
    m = vcat(Q.*M[1:r,:],M[r+1:end,:]);
    R = m \ b
    ea = abs(mean((Q-R')./Q));
	if ea < tol; return R; break; end
    Q = (Q+R')./2;
  end 
end

# ╔═╡ da63075a-0663-11eb-2557-272296908264
begin
	using Statistics, Printf
	# Valores de f, 𝜌 y D
	f = 0.005; 𝜌=1.23; D = 0.5;
	#    Longitudes de cada uno de los Tubos
	#    1 2 3 4 5 6 7 8 9 10
	L = [0 4 2 4 2 4 8 2 2  0 ];
	# Calulamos los valores de C correspondientes
	C = 16*f*L*𝜌/(2*(pi^2)*(D^5));
	# Planteamos la matriz de solucion
	#      Q2    Q3    Q4    Q5    Q6    Q7    Q8    Q9    Q10 
	N = [-C[2]  C[3]  C[4]    0     0     0     0   C[9]    0   ;
			0     0  -C[4]  C[5]  C[6]    0   C[8]    0     0   ;
			0     0     0     0  -C[6]  C[7]    0     0     0   ;
			1     1     0     0     0     0     0     0     0   ;
			0    -1     1     1     0     0     0     0     0   ;
			0     0     0    -1     1     1     0     0     0   ;
			0     0     0     0     1     1    -1     0     0   ;
			0     0     1     0     0     0     1    -1     0   ;
			1     0     0     0     0     0     0     1    -1  ;]
	# Vector de términos independientes
	b = [   0     0     0     1     0     0     0     0     0]';
	Q = Sol(N,b,0.1,0.01,3,100);
end

# ╔═╡ d25fa93e-0666-11eb-0ab3-ad9497f19aa8
md""" 
el valor inicial $Q_0$ es el mismo para todos los caudales $Q_i$ y se supuso inicialmente como $Q_i^* = 0.1 m^3/s$

Resolviendo el sistema los valores de los caudales $Q_i$  con $i=1,2,\cdots,10$ son:
"""

# ╔═╡ Cell order:
# ╟─7dbdf4d8-0587-11eb-0ce5-d397912d3eee
# ╟─119fe18e-0588-11eb-088b-8120e27e9433
# ╟─648d66d2-0588-11eb-2b95-4517fc61e165
# ╠═20e7ed8e-0593-11eb-1afb-7d77519ef9ee
# ╟─569072c0-0594-11eb-05e2-1985663cc0cc
# ╟─0c4f5008-0596-11eb-1db3-c14571aee5bb
# ╟─381b8ac4-0598-11eb-23bc-03663f91c2f3
# ╠═aa3d2694-0598-11eb-3b34-11a65f96e210
# ╠═76a9ce70-059a-11eb-13b5-13c94ec1c3e0
# ╟─cc74b18a-059a-11eb-3167-8525aeee8f3f
# ╟─684870e0-059d-11eb-3df9-21dba1e9dce7
# ╟─bbe25158-05ac-11eb-3cb6-ed5e87450d02
# ╟─583ec42c-05ae-11eb-1a75-051d43500d23
# ╠═095a0834-05af-11eb-2a36-d1d1e179e27d
# ╟─f25fb776-05b0-11eb-37bd-09b89dc26a2f
# ╠═5291f79e-05b1-11eb-0233-05ca4a69411b
# ╠═4c4d09fe-05b2-11eb-1609-b1e4adb4b50e
# ╟─a9445d72-05b4-11eb-0c50-e7b437caf6f2
# ╟─0ab3162a-05b5-11eb-0cc6-1bce617dfa96
# ╟─29365c30-05b6-11eb-331c-19a485128ce2
# ╟─6c699394-065a-11eb-26d7-1b895b06b2d4
# ╟─27dd8786-065e-11eb-0f4b-1ddcb3e6c6f9
# ╟─0fef8ab6-065e-11eb-3903-cfbc89172dad
# ╟─876e567e-0660-11eb-2906-91173842bfad
# ╠═aa576bf0-0663-11eb-1b89-35d73b118f51
# ╟─d25fa93e-0666-11eb-0ab3-ad9497f19aa8
# ╠═da63075a-0663-11eb-2557-272296908264
