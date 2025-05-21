This Python project simulates the motion of a body orbiting the Sun using two numerical methods: Forward Euler and Euler-Cromer.
It models a two-body gravitational system in three dimensions, likely representing a planet or asteroid orbiting the Sun, and generates plots to visualize the results.

Purpose:
The code solves the two-body problem, calculating the trajectory of an object under the gravitational influence of the Sun. It uses initial position and velocity data 
(likely for a celestial body like a planet or asteroid) and applies numerical integration to compute the object's path over time. The simulation compares two integration
methods—Forward Euler and Euler-Cromer—over 1-year and 3-year periods, with different time step sizes (dt), and visualizes the results.

Units and Constants:
Defines conversion factors for units like astronomical units (AU), meters, seconds, days, and years.
Sets universal constants: gravitational constant (G = 6.674e-11 m³ kg⁻¹ s⁻²) and the Sun's mass (m = 1.989e30 kg).

Initial Conditions:
Specifies initial positions (x0, y0, z0) in meters, converted from AU.
Specifies initial velocities (vx0, vy0, vz0) in meters per second, converted from AU per day.
The values suggest an orbit similar to a planet or asteroid, with a small z-coordinate indicating a near-planar orbit.

Time Parameters:
Simulates two time spans: 1 year and 3 years (converted to seconds).
Uses multiple time steps (dt = [0.1, 0.01, 0.001, 0.00001] years) to test numerical accuracy.

Numerical Methods:
Forward Euler Method (fc_solve_2body_fe):
A simple first-order numerical integration method.
Updates velocity and position iteratively based on gravitational acceleration.
Computes positions (x, y, z) and velocities (vx, vy, vz) over time.

Euler-Cromer Method (fc_solve_2body_ec):
A symplectic variant of the Euler method, improving energy conservation.
Updates velocity first, then uses the updated velocity to compute the new position.
Also computes positions and velocities in 3D.

Gravitational Model:
Models the Sun as a fixed central mass at the origin.
Uses Newton's law of gravitation to calculate acceleration: a = -GM * r / |r|^3, where r is the position vector.

Plotting:
Generates six plots for each method, time span, and time step:
Position vs. time (x vs t, y vs t, z vs t).
Position vs. position (z vs x, z vs y, y vs x), likely to show the orbital path.
Plots are created using matplotlib, with labels indicating the method, time step, and duration (1 or 3 years).

Loop Structure:
Iterates over each time step (dt) and computes trajectories for both methods over 1 and 3 years.
Calls plotting functions to visualize results for each case.

Outputs:
The project produces a series of plots showing:
How the x, y, and z coordinates evolve over time.
The orbital path in 2D projections (e.g., y vs x for the orbital plane).
Comparisons of the Forward Euler and Euler-Cromer methods for different time steps and durations.
The plots help assess the accuracy and stability of each method (e.g., smaller dt values yield more accurate results, and Euler-Cromer typically conserves energy better).

Likely Application for Further Study:
The initial conditions (e.g., position ~1 AU, velocities in AU/day) suggest the simulation models a body in a near-circular orbit around the Sun, possibly resembling Earth's orbit or an asteroid's.
The project likely serves as an educational tool, demonstrating numerical integration techniques in orbital mechanics, comparing the accuracy of Forward Euler (less stable) vs. Euler-Cromer (more stable for oscillatory systems like orbits).
