#define _USE_MATH_DEFINES
#include <iostream>
#include<stdlib.h>
#include<vector>
#include<cmath>
#include<fstream>

using namespace std;

//Equations in terms of double derivative of cartesian coordinates of the mass in terms of coordinates of the mass and centre of mass of the two bodies....
//CoM moving with constant velocity which may or may not be zero...
//Note to self: Remember to delete content of output text files when reducing n_iter...
//Mass values for Deimos-Mars system, used because of short orbital period,quick simulation, debugging. Very poor matching with measured values though, possibly because Pphobos not included.

const long double G = 6.67E-11;
const long double M = 6.41E23;
const long double m = 1.47E15;
const long double n_iter = 1000000;
const long double global_time_step = 0.1;



vector<vector<long double>>integrator(vector<long double>curr_positions, vector<long double>curr_velocities, long double time_step, long double time_sim,const long double eqn_const,vector<long double>curr_cm_position,vector<long double>curr_cm_velocity) {
    long double time = 0;

    vector<vector<long double>>positions;
    vector<long double>locations(3);
    locations = curr_positions;
    positions.push_back(locations);

    

    long double x_dot_total=0;
    long double y_dot_total=0;
    long double z_dot_total=0;

    long double x_dot_1,x_dot_2,x_dot_3,x_dot_4 = 0;
    long double y_dot_1,y_dot_2,y_dot_3,y_dot_4 = 0;
    long double z_dot_1,z_dot_2,z_dot_3,z_dot_4 = 0;

    while (1 == 1) {
        x_dot_total = 0;
        y_dot_total = 0;
        z_dot_total = 0;


        if (time > time_sim) {
            break;
        }

        //first derivative-LHS of time step

        x_dot_1 = curr_velocities[0];
        y_dot_1 = curr_velocities[1];
        z_dot_1 = curr_velocities[2];

        //second derivative-middle of time step

        x_dot_2 = x_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[0] - curr_cm_position[0]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));
        y_dot_2 = y_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[1] - curr_cm_position[1]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));
        z_dot_2 = z_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[2] - curr_cm_position[2]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));


        //third derivative-middle of time step

        x_dot_3 = x_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[0] + x_dot_1 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_1 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_1 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_1 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));
        y_dot_3 = y_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[1] + y_dot_1 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_1 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_1 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_1 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));
        z_dot_3 = z_dot_1 - 0.5 * time_step * eqn_const * (curr_positions[2] + z_dot_1 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_1 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_1 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_1 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));

        //fourth derivative-RHS of time step


        x_dot_4 = x_dot_1 - time_step * eqn_const * (curr_positions[0] + x_dot_2 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_2 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_2 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_2 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));
        y_dot_4 = y_dot_1 - time_step * eqn_const * (curr_positions[1] + y_dot_2 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_2 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_2 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_2 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));
        z_dot_4 = z_dot_1 - time_step * eqn_const * (curr_positions[2] + z_dot_2 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[0] * time_step / 2) / (pow(pow((curr_positions[0] + x_dot_2 * time_step / 2 - curr_cm_position[0] - curr_cm_velocity[0] * time_step / 2), 2) + pow((curr_positions[1] + y_dot_2 * time_step / 2 - curr_cm_position[1] - curr_cm_velocity[1] * time_step / 2), 2) + pow((curr_positions[2] + z_dot_2 * time_step / 2 - curr_cm_position[2] - curr_cm_velocity[2] * time_step / 2), 2), 1.5));

        //Runge-kutta method

        x_dot_total = (x_dot_1 + 2 * x_dot_2 + 2 * x_dot_3 + x_dot_4) / 6;
        y_dot_total = (y_dot_1 + 2 * y_dot_2 + 2 * y_dot_3 + y_dot_4) / 6;
        z_dot_total = (z_dot_1 + 2 * z_dot_2 + 2 * z_dot_3 + z_dot_4) / 6;

        curr_positions[0] = curr_positions[0] + x_dot_total * time_step;
        curr_positions[1] = curr_positions[1] + y_dot_total * time_step;
        curr_positions[2] = curr_positions[2] + z_dot_total * time_step;

        locations = curr_positions;

        positions.push_back(locations);

        time = time + time_step;

        curr_cm_position[0] = curr_cm_position[0] + curr_cm_velocity[0] * time_step;
        curr_cm_position[1] = curr_cm_position[1] + curr_cm_velocity[1] * time_step;
        curr_cm_position[2] = curr_cm_position[2] + curr_cm_velocity[2] * time_step;

        curr_velocities[0] = curr_velocities[0] - time_step * eqn_const * (curr_positions[0] - curr_cm_position[0]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));
        curr_velocities[1] = curr_velocities[1] - time_step * eqn_const * (curr_positions[1] - curr_cm_position[1]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));
        curr_velocities[2] = curr_velocities[2] - time_step * eqn_const * (curr_positions[2] - curr_cm_position[2]) / (pow(pow((curr_positions[0] - curr_cm_position[0]), 2) + pow((curr_positions[1] - curr_cm_position[1]), 2) + pow((curr_positions[2] - curr_cm_position[2]), 2), 1.5));


    }

    return positions;
   
}

long double vector_mod(vector<long double>A) {

    return pow(pow(A[0], 2) + pow(A[1], 2) + pow(A[2], 2), 0.5);

}

vector<long double>cross_product(vector<long double>A, vector<long double>B) {
    vector<long double>C(3);
    C[0]= A[1] * B[2] - A[2] * B[1];
    C[1]= A[2] * B[0] - A[0] * B[2];
    C[2]= A[0] * B[1] - A[1] * B[0];

    return C;
}

long double dot_product(vector<long double>A, vector<long double>B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];

}

vector<vector<long double>> mat_mul_3_3(vector<vector<long double>>A, vector<vector<long double>>B) {
    vector <vector<long double>>C;
    vector<long double>c_row;
    long double c_element;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c_element = 0;
            for (int k = 0; k < 3; k++) {
                c_element = c_element + A[i][k] * B[k][j];
            }
            c_row.push_back(c_element);
        }
        C.push_back(c_row);
        c_row.clear();
    }


    return C;
}

vector<long double>mat_mul_3_1(vector<vector<long double>>A, vector<long double>B) {
    vector<long double>C;
    long double temp_val = 0;
    for (int i = 0; i < 3; i++) {
        temp_val = 0;
        for (int j = 0; j < 3; j++) {
            temp_val = temp_val + A[i][j] * B[j];

        }
        C.push_back(temp_val);
        
    }

    return C;

}

long double det(vector<vector<long double>>A) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) + A[0][1] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

}

vector<vector<long double>>make_t(vector<vector<long double>>A) {
    vector<vector<long double>>B;
    vector<long double> B_row(3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            B_row[j] = A[j][i];
        }
        B.push_back(B_row);
    }

    return B;

}

vector<vector<long double>>scalar_multiple(vector<vector<long double>>A,long double mul) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            A[i][j] = A[i][j] * mul;
        }
    }
    return A;
}

void testing_module(vector<long double>curr_positions_1, vector<long double>curr_velocities_1, vector<long double>curr_positions_2,vector<long double>curr_velocities_2,const long double eqn_const_1,const long double eqn_const_2, vector<long double>curr_cm_position, long double time_sim, vector<long double>curr_cm_velocity) {
    //Analytical solution found and sent to text files to evaluate numerical solution
    vector<long double>temp;
    vector<long double>temp1;
    vector<long double>temp2;
    long double temp_val = 0;


    //All values in COM frame
    vector<long double>relative_pos_1(3);
    relative_pos_1[0] = curr_positions_1[0] - curr_cm_position[0];
    relative_pos_1[1] = curr_positions_1[1] - curr_cm_position[1];
    relative_pos_1[2] = curr_positions_1[2] - curr_cm_position[2];

    vector<long double>relative_pos_2(3);
    relative_pos_2[0] = curr_positions_2[0] - curr_cm_position[0];
    relative_pos_2[1] = curr_positions_2[1] - curr_cm_position[1];
    relative_pos_2[2] = curr_positions_2[2] - curr_cm_position[2];


    vector<long double>relative_vels_1(3);
    relative_vels_1[0] = curr_velocities_1[0] - curr_cm_velocity[0];
    relative_vels_1[1] = curr_velocities_1[1] - curr_cm_velocity[1];
    relative_vels_1[2] = curr_velocities_1[2] - curr_cm_velocity[2];

    vector<long double>relative_vels_2(3);
    relative_vels_2[0] = curr_velocities_2[0] - curr_cm_velocity[0];
    relative_vels_2[1] = curr_velocities_2[1] - curr_cm_velocity[1];
    relative_vels_2[2] = curr_velocities_2[2] - curr_cm_velocity[2];


    //Calculating orbital parameters
    long double energy = 0.5 * (pow(vector_mod(relative_vels_1), 2)) - eqn_const_1 / (vector_mod(relative_pos_1));
    long double h = vector_mod(cross_product(relative_pos_1, relative_vels_1));
    long double e = sqrt(1+2*h*h*energy/(eqn_const_1*eqn_const_1));
    long double a = h * h / (eqn_const_1 * (1 - e * e));
    long double theta = acos((h*h/(eqn_const_1*vector_mod(relative_pos_1)) - 1) / e);
    long double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
    long double time = a * sqrt(a / eqn_const_1) * (E - e*sin(E));
    long double time_0 = time;
    long double theta_increment = 0.01;
    long double r_mod;
    vector<vector<long double>>positions;
    positions.push_back(curr_positions_1);
    vector<long double>coordinates(3);

    //Angles made by orbital plane, used to convert values produced by solution to ground frame.
    vector<long double>orbit_plane_normal;
    vector<vector<long double>>I;

    orbit_plane_normal = cross_product(relative_vels_1,relative_vels_2);
    
    long double mod_orbit_plane_normal = vector_mod(orbit_plane_normal);
    orbit_plane_normal[0] = orbit_plane_normal[0] / mod_orbit_plane_normal;
    orbit_plane_normal[1] = orbit_plane_normal[1] / mod_orbit_plane_normal;
    orbit_plane_normal[2] = orbit_plane_normal[2] / mod_orbit_plane_normal;

    long double plane_angle_x, plane_angle_y, plane_angle_z;
    
    temp.push_back(1);
    temp.push_back(0);
    temp.push_back(0);
    I.push_back(temp);
    plane_angle_x = M_PI_2 - (dot_product(orbit_plane_normal, temp) / (mod_orbit_plane_normal));
    temp[0] = 0;
    temp[1] = 1;
    I.push_back(temp);
    plane_angle_y = M_PI_2 - (dot_product(orbit_plane_normal, temp) / (mod_orbit_plane_normal));
    temp[1] = 0;
    temp[2] = 1;
    I.push_back(temp);
    plane_angle_z = M_PI_2 - (dot_product(orbit_plane_normal, temp) / (mod_orbit_plane_normal));

    long double first_rot = M_PI_2 - plane_angle_z;
    long double second_rot;
    
   
    while(theta<2*M_PI ){
        r_mod = h * h / (eqn_const_1 * (1 + e * cos(theta)));
        coordinates[0] = r_mod * cos(theta);
        coordinates[1] = r_mod * sin(theta);
        coordinates[2] = 0;

        //transforming coordinates by aligning z axis first

        vector<long double>rot_axis = cross_product(orbit_plane_normal, temp);
        long double mod_rot_axis = vector_mod(rot_axis);
        rot_axis[0] = rot_axis[0] / mod_rot_axis;
        rot_axis[1] = rot_axis[1] / mod_rot_axis;
        rot_axis[2] = rot_axis[2] / mod_rot_axis;

        vector<vector<long double>>rot_axis_into_t;
        temp.clear();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                temp.push_back(rot_axis[i] * rot_axis[j]);
                
            }
            rot_axis_into_t.push_back(temp);
            temp.clear();
        }


        temp.clear();
        temp = mat_mul_3_1(scalar_multiple(I, cos(first_rot)), coordinates);
        temp1 = mat_mul_3_1(scalar_multiple(rot_axis_into_t, (1 - cos(first_rot))),coordinates);
        temp2 = cross_product(rot_axis, coordinates);
        temp2[0] = temp2[0] * sin(first_rot);
        temp2[1] = temp2[1] * sin(first_rot);
        temp2[2] = temp2[2] * sin(first_rot);
        
        temp[0] = temp[0] + temp1[0] + temp2[0];
        temp[1] = temp[1] + temp1[1] + temp2[1];
        temp[2] = temp[2] + temp1[2] + temp2[2];

        coordinates = temp;
        
        temp.clear();
        temp1.clear();
        temp2.clear();

        temp.push_back(1);
        temp.push_back(0);
        temp.push_back(0);

        second_rot = dot_product(coordinates, temp) / vector_mod(coordinates);

        temp = coordinates;

        coordinates[0] = temp[0] * cos(second_rot) - temp[1] * sin(second_rot);
        coordinates[1] = temp[0] * sin(second_rot) + temp[1] * cos(second_rot);
        
        temp.clear();

        coordinates[0] = coordinates[0] + curr_cm_position[0] + curr_cm_velocity[0] * time;
        coordinates[1] = coordinates[1] + curr_cm_position[1] + curr_cm_velocity[1] * time;
        coordinates[2] = coordinates[2] + curr_cm_position[2] + curr_cm_velocity[2] * time;

        positions.push_back(coordinates);

        theta = theta + theta_increment;
        E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(theta / 2));
        time = a * sqrt(a / eqn_const_1) * (E - e * sin(E));

    }

    ofstream testdatafile1("test1.txt");
    for (int i = 0; i < positions.size(); i++) {
        testdatafile1 << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;

    }
    testdatafile1.close();

}






int main(void) {
//Subscript 1 for m, 2 for M. In input, m is taken first.
    vector<long double>input_positions_1(3);
    vector<long double>input_velocities_1(3);
    vector<long double>input_positions_2(3);
    vector<long double>input_velocities_2(3);
    vector<long double>cm_pos(3);
    vector<long double>cm_vel(3);
    long double input_number;
    for (int i = 0; i < 3; i++) {
        cin >> input_number;
        input_positions_1[i]=input_number;
        
    }

    for (int i = 0; i < 3; i++) {
        cin >> input_number;
        input_velocities_1[i]=input_number;

    }

    for (int i = 0; i < 3; i++) {
        cin >> input_number;
        input_positions_2[i] = input_number;

    }

    for (int i = 0; i < 3; i++) {
        cin >> input_number;
        input_velocities_2[i] = input_number;

    }

    for (int i = 0; i < 3; i++) {
        cm_pos[i] = (m * input_positions_1[i] + M * input_positions_2[i]) / (M + m);
    }

    for (int i = 0; i < 3; i++) {
        cm_vel[i] = (m * input_velocities_1[i] + M * input_velocities_2[i]) / (M + m);
    }
 

    vector<vector<long double>>data_1 = integrator(input_positions_1,input_velocities_1,global_time_step, n_iter,G*M*M*M/((M+m)*(M+m)),cm_pos,cm_vel);
    vector<vector<long double>>data_2 = integrator(input_positions_2, input_velocities_2, global_time_step, n_iter, G * m * m * m / ((M + m) * (M + m)), cm_pos, cm_vel);

    //Writing output coordinates in x-y-z to text files for plotting with gnuplot.Plotting was done via terminal in the relevant folder.

    ofstream datafile1("orbitvalues1.txt");
    for (int i = 0; i < n_iter; i++) {
        datafile1 << data_1[i][0] << " " << data_1[i][1] << " " << data_1[i][2]<<endl;

    }
    datafile1.close();

    ofstream datafile2("orbitvalues2.txt");
    for (int i = 0; i < n_iter; i++) {
        datafile2 << data_2[i][0] << " " << data_2[i][1] << " " << data_2[i][2] << endl;

    }
    datafile2.close();

    //Test module, use if orbit is inclined, neither of the masses has initial velocity zero.
    //testing_module(input_positions_1,input_velocities_1,input_positions_2,input_velocities_2, G * M * M * M / ((M + m) * (M + m)), G * m * m * m / ((M + m) * (M + m)),cm_pos, n_iter, cm_vel);



    return 0;
}