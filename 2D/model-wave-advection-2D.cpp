#include <iostream>
#include <fstream>
#include <cmath>

int main() {

  // Define the parameters
  double u = 0.5; // speed of the pulse in x direction
  double v = 0.5; // speed of the pulse in y direction
  double sigma_x = 0.05; // width of the pulse in x direction
  double sigma_y = 0.05; // width of the pulse in y direction
  double xc = 0.5; // center location of the pulse in x direction
  double yc = 0.5; // center location of the pulse in y direction
  double dx = 0.01; // spatial x step size
  double dy = 0.01; // spatial y step size
  double dt = 0.01; // time step size
  int nx = 100; // number of spatial grid points in x
  int ny = 100; // number of spatial grid points in y
  int nt = 200; // number of time steps
  
  // Define variables
  double x[nx]; // array for x grid
  double y[ny]; // array for y grid
  double t[nt]; // array for time
  double psi[nx][ny]; // array for wave amplitude
  double psi_new[nx][ny]; // array for previous time step of psi
  double x_diff; // offset of pulse from the center in y
  double y_diff; // offset of pulse from the center in x
  double exponent_x; // exponent of pulse in x
  double exponent_y; // exponent of pulse in y
  
  // Open the files for writing
  std::ofstream psi_file("data/psi.dat");
  std::ofstream x_file("data/x.dat");
  std::ofstream y_file("data/y.dat");
  std::ofstream t_file("data/t.dat");

  // Set up the x array
  for (int i=0; i<nx; i++) {

    // Define the array value
    x[i] = i*dx;

    // Write to file
    x_file << x[i] << " ";
    
  }
  x_file << std::endl;
  x_file.close();
    
  // Set up the y array
  for (int j=0; j<ny; j++) {

    // Define the array value
    y[j] = j*dy;

    // Write to file
    y_file << y[j] << " ";
    
  }
  y_file << std::endl;
  y_file.close();

  // Set up the t array
  for (int n=0; n<nt; n++) {

    // Define the array value
    t[n] = n*dt;

    // Write to file
    t_file << t[n] << " ";
    
  }
  t_file << std::endl;
  t_file.close();
    
  // Loop through x and y and initialize the pulse at t=0
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
    
      // Get the pulse offset from 0 for x and y
      x_diff = x[i]-xc;
      y_diff = y[j]-yc;

      // Build the x and y parts of the exponent
      exponent_x = -(x_diff*x_diff)/(2.0*sigma_x*sigma_x);
      exponent_y = -(y_diff*y_diff)/(2.0*sigma_y*sigma_y);

      // Build the Gaussian pulse
      psi[i][j] = exp(exponent_x + exponent_y);
     
    }
  }

  // Apply the y boundary conditions
  for (int i=0; i<nx; i++) {
    psi[i][0] = psi[i][ny-1];
  }

  // Apply the x boundary conditions
  for (int j=0; j<ny; j++) {
    psi[0][j] = psi[nx-1][j];
  }

  // Write the initial psi to file
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      psi_file << psi[i][j] << " ";
    }
  }
  psi_file << std::endl;
  
  // Time evolution of the pulse
  for (int n=1; n<nt; n++) {
    
    // Update psi_new
    for (int i=1; i<nx; i++) {
      for (int j=1; j<ny; j++) {
	psi_new[i][j] = psi[i][j] - u*dt/dx*(psi[i][j]-psi[i-1][j]) - v*dt/dy*(psi[i][j]-psi[i][j-1]);
      }
    }

    // Apply the y boundary conditions
    for (int i=0; i<nx; i++) {
      psi_new[i][0] = psi_new[i][ny-1];
    }

    // Apply the x boundary conditions
    for (int j=0; j<ny; j++) {
      psi_new[0][j] = psi_new[nx-1][j];
    }

    // Update psi and write to file
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {

	// Update psi
	psi[i][j] = psi_new[i][j];

	// Write psi to output
	psi_file << psi[i][j] << " ";
	
      }
    }
    psi_file << std::endl;
   
  }
  // close the file
  psi_file.close();
  
  return 0;
  
}
