#include <iostream>
#include <fstream>
#include <cmath>

// Define a function to return the Gaussian pulse
double gaussian_pulse(double x, double xc, double sigma) {
    double exponent = -((x - xc) * (x - xc)) / (2.0 * sigma * sigma);
    return exp(exponent);
}

int main() {

  // Define the parameters
  double c = 1.0; // speed of the pulse
  double sigma = 0.1; // width of the pulse
  double xc = 0.5; // center location of the pulse
  double dx = 0.01; // spatial step size
  double dt = 0.01; // time step size
  int nx = 100; // number of spatial grid points
  int nt = 100; // number of time steps
  
  // Define variables
  double x, t; // space and time grid values
  double u[nx]; // array for wave amplitude
  double uold[nx]; // array for previous time step of u
  
  // Open the files for writing
  std::ofstream ufile("u.dat");
  std::ofstream xfile("x.dat");
  std::ofstream tfile("t.dat");
  
  // Initialize the pulse at t=0
  for (int i = 0; i < nx; i++) {

    // Get current x value and write to file
    x = i * dx;
    xfile << x << " ";

    // Generate the Gaussian pulse and write to file
    u[i] = gaussian_pulse(x, xc, sigma);
    ufile << u[i] << " ";

  }
  xfile << std::endl;
  xfile.close();
  ufile << std::endl;
    
  // Apply the boundary conditions
  u[0] = u[nx-1];
  
  // Time evolution of the pulse
  for (int n = 1; n < nt; n++) {
    
    // Get the current t value and write to file
    t = n * dt; 
    tfile << t << " ";

    // Set up uold
    for (int i = 0; i < nx; i++) {
      uold[i] = u[i];
    }
    
    // Update u
    for (int i = 1; i < nx; i++) {
      u[i] = uold[i] - c * dt / dx * (uold[i] - uold[i-1]);
    }

    // Apply the boundary conditions
    u[0] = u[nx-1];
    
    // Write the pulse to the file at the current time
    for (int i = 0; i < nx; i++) {
      ufile << u[i] << " ";
    }
    ufile << std::endl;
    
  }
  tfile << std::endl;
  tfile.close();
  
  // close the file
  ufile.close();
  
  return 0;
  
}
