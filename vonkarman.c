// Iniciacion de la simulacion de von Karman, importando las librerias necesarias
#include "embed.h" //Librería para el objeto en el que incide el fluido
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
const double Reynolds = 160.;
int maxlevel = 9;
face vector muv[];


//Dominio del problema, 8 unidades de largo y centrado en el origen:
int main()
{
  L0 = 8. [1];
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  //Manejo del la precision de la malla de refinamiento
  display_control (Reynolds, 10, 1000);
  display_control (maxlevel, 6, 12);
  run(); 
}


double D = 0.125, U0 = 1.; //Diametro del cilindro y velocidad de entrada, respectivamente
//Seteo una viscosidad constante en base a los datos establecidos para el problema
event properties (i++)
{
  foreach_face()
  muv.x[] = fm.x[]*D*U0/Reynolds;
}


//Condiciones de frontera del fluido:
u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);


//Condiciones de frontera del cilindro:
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
  solid (cs, fs, sqrt(sq(x) + sq(y)) - D/2.); //Diametro del cilindro de 0.125
  foreach()
  u.x[] = cs[] ? U0 : 0.; //Seteo el inicio del campo de velocidad
}


//Checkeo el número de iteraciones
event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);


//Produzco algunas animaciones iniciales de la simulacion:
event movies (i += 4; t <= 15.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      linear = false, min = 0, max = 1, mask = m);
}


//Adaptaciones en base a los posibles errores del modelamiento y sus componentes:
event adapt (i++) {
  adapt_wavelet ({cs,u,f}, {1e-2,3e-2,3e-2,3e-2}, maxlevel, 4);
}