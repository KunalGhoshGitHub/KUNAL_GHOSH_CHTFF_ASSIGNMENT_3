#include <iostream>
#include <cassert>
#include<algorithm>
#include<vector>
#include<fstream>
#include<cmath>
#include<string>
#include<unistd.h>
#include<functional>

using namespace std;

typedef double real;

double const pi = 4.0*atan(1.0);

// Generate 1D grid
void generate_1D_grid(vector<real> &x, real X0, real Xl, real N)
{
    for (int i = 0; X0 + ((i/N)*(Xl-X0)) <= Xl;i++)
    {
        x.push_back(X0 + ((i/N)*(Xl-X0)));
        //cout<<X0 + ((i/N)*(Xl-X0))<<endl;
    }
}

// This function calculates the transpose of a given matrix
void matrix_transpose(vector<vector<real>>& A, vector<vector<double>> &A_Transpose)
{
    int rows = A.size();
    int cols = A[0].size();

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            A_Transpose[col][row] = A[row][col];
        }
    }
}

// This function calculates the eigen values , eigen vector matrix and inverse of eigen vector matrix
void compute_eigen_value_vector_inverse(vector<vector<real>> &P, vector<vector<real>> &PINV, vector<real> &eigen_values, vector<real> &x,real Nx,real hx,real nu)
{
    // Set up eigenvalues and matrices
    int i, j;
	for(i = 1 ; i < Nx ; i++) {
		eigen_values.push_back(0.5-(-nu*4.0*sin(pi*0.5*x[i])*sin(pi*0.5*x[i])));
		for(j = 1 ; j < Nx ; j++) {
			P[i-1][j-1] = sin(i*pi*x[j])/sqrt((Nx-1.0)/2.0) ;
			PINV[i-1][j-1] = sin(j*pi*x[i])/sqrt((Nx-1.0)/2.0) ;
		}
	}
}

// This function calculate the F of the equation at time t
void compute_F(vector<vector<real>> &F, vector<real> &x, vector<real> &y, real h, real t)
{
    int rows = F.size();
    int cols = F[0].size();
    real pi = 4.0*atan(1);

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            F[row][col] = h*h*sin(2.0*pi*x[row+1])*sin(2.0*pi*y[col+1])*sin(2.0*pi*t);
        }
    }
}

// This function calculate the analytical solution
void compute_analytical_sol(vector<vector<real>> &sol, vector<real> &x, vector<real> &y, real h, real t)
{
    int rows = sol.size();
    int cols = sol[0].size();
    real pi = 4.0*atan(1);
    real temp = (1.0 + (16*pi*pi));

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            sol[row][col] = sin(2.0*pi*x[row+1])*sin(2.0*pi*y[col+1])*((exp(-8.0*pi*pi*t)/(2*pi*temp)) + (2.0*sin(2*pi*t)/temp) - (cos(2.0*pi*t)/(2*pi*temp)));
        }
    //cout<<x[row+1]<<endl;
    }
}



// This function computes the product of two given matrix
void matrix_multiply(vector<vector<real>> &A, vector<vector<double>> &B, vector<vector<double>> &C)
{
    int m = A.size();
    int n = B.size();
    int p = B[0].size();

    // vector<vector<double>> C(m, vector<double>(p, 0));

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    //return C;
}

// This function computes the addition of two given matrix
void matrix_addition(vector<vector<real>> &A, vector<vector<double>> &B, vector<vector<double>> &C)
{
    int m = A.size();
    int n = A[0].size();
    // int n = B.size();
    // int p = B[0].size();

    // vector<vector<double>> C(m, vector<double>(p, 0));

    for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                C[i][j] = A[i][j] + B[i][j];
            }
        }

    //return C;
}

// This function computes the addition of two given matrix
void matrix_value_transfer(vector<vector<real>> &A, vector<vector<double>> &B)
{
    int m = A.size();
    int n = A[0].size();
    // int n = B.size();
    // int p = B[0].size();

    // vector<vector<double>> C(m, vector<double>(p, 0));

    for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                A[i][j] = B[i][j];
            }
        }

    //return C;
}

// This function computes the scalar multiplication of a given matrix
void matrix_scalar_multiply(vector<vector<real>> &A,real nu, vector<vector<double>> &C)
{
    int m = A.size();
    int n = A[0].size();

    for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                C[i][j] = nu*A[i][j];
            }
        }

    //return C;
}

// This function computes the V matrix
void compute_V(vector<vector <real>> &F_tilde, vector<real> &eigen_values, vector<vector<real>> &V)
{
    int rows = F_tilde.size();
    int cols = F_tilde[0].size();

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            V[row][col] = (F_tilde[row][col])/(eigen_values[row] + eigen_values[col]);
        }
    }
}

// This function computes the V matrix
void compute_D(vector<vector <real>> &D)
{
    int rows = D.size();
    int cols = D[0].size();

    for (int row = 0; row < rows; row++)
    {
            D[row][row] = -2.0;
    }
    for (int row = 0; row < rows-1; row++)
    {
            D[row][row+1] = 1;
    }
    for (int row = 1; row < rows; row++)
    {
            D[row][row-1] = 1;
    }
}

// Function to save a vector<double> to a file of given name
void write_to_file(vector<real> &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    cout<<"Output file "<< str <<" is opened."<<endl;

    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<scientific<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// Function to save a vector<string> to a file of given name
void write_to_file(vector<string> &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    cout<<"Output file "<< str <<" is opened."<<endl;

    // Writing the vector values to the file in scientific notation
    for(int i=0;i<u.size();i++)
    {
        file<<u[i]<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

// Function to save a vector<vector <double>> (MATRIX) to a file of given name
void write_to_file(vector<vector <real> > &u, string str)
{
    ofstream file;
    // Open the file
    file.open(str);

    // If the file is NOT opened
    if( !file )
    {
        // Error Message if the file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }

    cout<<"Output file "<< str <<" is opened."<<endl;

    int rows = u.size();
    int cols = u[0].size();

    // Writing the vector values to the file in scientific notation

    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            file<<u[row][col]<<scientific<<",";
        }
        file<<endl;
    }

    // Closing the file
    file.close();

    cout<<"Output file "<< str <<" is closed."<<endl;
    cout<<endl;
}

double calculate_average(const std::vector<real>& v) {
    real sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    double avg = static_cast<double>(sum) / v.size();
    return avg;
}

int calculate_max(const std::vector<real>& v) {
    real max = v[0];
    for (int i = 1; i < v.size(); i++) {
        if (v[i] > max) {
            max = v[i];
        }
    }
    return max;
}


int main()
{
    cout << "Hello world!" << endl;
    real dt,nu;

    real X0,Xl,Y0,Yl,h;
    int N=100;

    // Given the domain of the solution
    X0 = 0;
    Xl = 1;
    Y0 = 0;
    Yl = 1;

    h = 1.0/N;

    dt = 0.005;
    nu = dt/(2.0*h*h);

    nu = dt/(2.0*h*h);
    int rows, cols;
    rows = N-1;
    cols = N-1;


    vector<real> x;
    vector<real> y;
    vector<vector <real>> U_extended(N,vector<real> (N,0.0));

    // Generating the Grids
    generate_1D_grid(x,X0,Xl,N);
    generate_1D_grid(y,Y0,Yl,N);

    vector<vector <real> > P(rows,vector<real> (cols,0.0));
    vector<vector <real> > P_INV(rows,vector<real> (cols,0.0));
    vector<real> eigen_values;
    // Computing the eigen values, eigen vectors and the inverse of the matrix of eigen vectors
    compute_eigen_value_vector_inverse(P,P_INV,eigen_values,x,N,h,nu);

    write_to_file(eigen_values,"Eigen_Value_.csv");

    // Vector to store the names of all the output files generated
    vector<string> Output_file_names;
    vector<string> Analytical_Output_file_names;

    vector<vector <real> > U_n(rows,vector<real> (cols,0.0));
    vector<vector <real> > Analytical_U_n(rows,vector<real> (cols,0.0));

    vector<real> L2_error_Vector, max_error_vector;

    real t_end;
    t_end = 1.0;
    t_end = t_end+dt;
    for (real t = 0.0; t <= t_end; t = t+dt){
        cout<<"Time = "<<t<<endl;
        vector<vector <real> > F_n(rows,vector<real> (cols,0.0));
        compute_F(F_n,x,y,h,t);

        vector<vector <real> > F_n_plus_1(rows,vector<real> (cols,0.0));
        compute_F(F_n_plus_1,x,y,h,t+dt);

        vector<vector <real> > D(rows,vector<real> (cols,0.0));
        compute_D(D);

        vector<vector <real> > temp(rows,vector<real> (cols,0.0));
        matrix_multiply(U_n,D,temp);

        vector<vector <real> > temp_1(rows,vector<real> (cols,0.0));
        matrix_multiply(D,U_n,temp_1);

        vector<vector <real> > temp_2(rows,vector<real> (cols,0.0));
        matrix_addition(F_n,F_n_plus_1,temp_2);

        vector<vector <real> > temp_3(rows,vector<real> (cols,0.0));
        matrix_addition(temp_2,temp_1,temp_3);

        vector<vector <real> > temp_4(rows,vector<real> (cols,0.0));
        matrix_addition(temp_3,temp,temp_4);

        vector<vector <real> > temp_5(rows,vector<real> (cols,0.0));
        matrix_scalar_multiply(temp_4,nu,temp_5);

        vector<vector <real> > temp_6(rows,vector<real> (cols,0.0));
        matrix_addition(U_n,temp_5,temp_6);

        vector<vector <real> > temp_7(rows,vector<real> (cols,0.0));
        matrix_multiply(temp_6,P,temp_7);

        vector<vector <real> > F_tilde(rows,vector<real> (cols,0.0));
        matrix_multiply(P_INV,temp_7,F_tilde);


        vector<vector <real> > V(rows,vector<real> (cols,0.0));
        compute_V(F_tilde,eigen_values,V);

        vector<vector <real> > temp_9(rows,vector<real> (cols,0.0));
        matrix_multiply(P,V,temp_9);

        vector<vector <real> > U_n_plus_1(rows,vector<real> (cols,0.0));
        matrix_multiply(temp_9,P_INV,U_n_plus_1);

        compute_analytical_sol(Analytical_U_n,x,y,h,t);

        // Writing the Numerical Solution to the file
        write_to_file(U_n,"Q_4_Numerical_Solution_at_time_"+to_string(t)+"_.csv");

        // Writing the Analytical Solution to the file
        write_to_file(Analytical_U_n,"Q_4_Analytical_Solution_at_time_"+to_string(t)+"_.csv");

        Output_file_names.push_back("Q_4_Numerical_Solution_at_time_"+to_string(t)+"_.csv");

        Analytical_Output_file_names.push_back("Q_4_Analytical_Solution_at_time_"+to_string(t)+"_.csv");

        write_to_file(Output_file_names,"Output_file_names.csv");

        write_to_file(Analytical_Output_file_names,"Analytical_Output_file_names.csv");

        real Max_Error, L2_Error;
        Max_Error = L2_Error = 0.0 ;


        // Write a Output file showing values of x,y, Numerical Solution and Analytical Solution
        ofstream File("Q_4_Output_"+to_string(N)+"_.dat", ios::out) ;
        File.flags( ios::dec | ios::scientific );
        File.precision(16) ;
        if(!File) {cerr<< "Error: Output file couldnot be opened.\n";}

        File << "TITLE = Flow" << endl << "VARIABLES = X, Y, u, Exact " << endl;
        File << "Zone T = psi I = " << N+1 << " J = " << N+1 << endl ;
        int i,j;

        for(i = 0 ; i < N-1 ; i++) {
            for(j = 0 ; j < N-1 ; j++) {
                if( fabs(U_n[i][j] - Analytical_U_n[i][j] ) > Max_Error)
                {
                    Max_Error = fabs( U_n[i][j] - Analytical_U_n[i][j] );
                }
                L2_Error += (U_n[i][j] - Analytical_U_n[i][j])*(U_n[i][j] - Analytical_U_n[i][j] )/ ( ( N+1.0 )*(N+1.0) ) ;

                File << x[i] << "\t" << y[j] << "\t" << U_n[i][j] << "\t" << Analytical_U_n[i][j] << endl ;
            }
        }


        L2_Error = sqrt(L2_Error) ;
        File.close() ;
        // Printing the value of N
        cout<<"For N = "<<N<<endl;

        L2_error_Vector.push_back(L2_Error);
        max_error_vector.push_back(Max_Error);

        matrix_value_transfer(U_n,U_n_plus_1);

    }

    vector<real> L2_Error;
    vector<real> Max_Error;

    L2_Error.push_back(calculate_average(L2_error_Vector));

    Max_Error.push_back(calculate_max(max_error_vector));

    // Writing the L2 Errors in a file
    write_to_file(L2_Error,"Q_4_L2_Error_h_"+to_string(h)+"_dt_"+to_string(dt)+".csv");

    // Writing the Max Errors in a file
    write_to_file(Max_Error,"Q_4_max_error_h_"+to_string(h)+"_dt_"+to_string(dt)+".csv");

    return 0;
}
