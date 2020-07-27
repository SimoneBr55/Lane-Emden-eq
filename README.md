# lane-emden-eq

This is the script in Fortran, I used in the course 'Astronomical Computing' exam (translated from "Calcolo per l'Astronomia") I took in July 2020.

It basically 

  1. integrates the Lane-Emden Equation numerically by using a 4th Order Runge-Kutta algorithm;
  2. finds the $\xi_0$ (radius of the star) by a bisection algorithm and
  3. computes the critical mass of Chandrasekhar.
 
# Astrophysical Problem

Lane-Emden Equation (LEE) is a dimensionless form of the Poisson Equation which studies the gravitational potential of a politropic fluid, under the condition of it begin spherical symmetrical and self-gravitating.

Mathematically, LEE is a 2nd Order Ordinary Differential Equation (ODE2) and is useful in Astrophysics to model the internal structure of a star.

It is valid under two main assumptions:

  1. the gas forming the stellar structure must be in hydrostatic equilibrium with gravity;
  2. the pressure of the gas and its density must be related by a politropic equation of state.
  
If these conditions are met, we can write a dimensionless equation linking the density (or pressure) with the radius of the star (or the structure, in general).

1/&xi; d/d&xi; ( &xi; <sup>2</sup> d&theta; / d&xi; ) = - &theta; <sup>n</sup>

I have put a PDF document in this repo to explain more extensively the astrophyisical problem and the mathematic derivation of the LEE. It is now available only in Italian. ASAP, I will translate it in English.
