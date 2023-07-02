#include<iostream>
#include<stdlib.h>
#include<vector>
#include<fstream>
#include<cmath>


double find_fast_sqrt(double my_number) {
	double answer;
	long long manipulated_input = *(long long *) &my_number; //Put double in address of long long, both 64 bits
	manipulated_input = 0x5FE6F7CED916872B - (manipulated_input >> 1); //Appropriate scaling for double-11 bits exponent, 52 bits mantissa, 1 sign bit.
	answer = *(double*)&manipulated_input; //Back to double
	//Newton iteration, two times
	answer = (answer * 3 - my_number * answer * answer * answer) / 2;
	answer = (answer * 3 - my_number * answer * answer * answer) / 2;
	return answer;

}

class celestial_body { //n instances of this class will be used, each containing position, momentum data and history
	double mass;
	std::vector<double>pos;
	std::vector<double>mom;
	std::vector<double>provisional_pos;
	std::vector<double>provisional_mom;
	std::vector<std::vector<double>>pos_history;
	std::vector<std::vector<double>>mom_history; 

	friend class n_body_data;
	friend class n_body_solution;
	friend class derivative_x;
	friend class derivative_y;
	friend class derivative_z;
	friend class derivative_px;
	friend class derivative_py;
	friend class derivative_pz;

};

class n_body_data {            //Contains the vector of all celestial_body instances as well as time, time_step, G, etc
	friend class n_body_solution;  //Calls n_body_data and executes the simulation
	friend class derivative_x;
	friend class derivative_y;
	friend class derivative_z;
	friend class derivative_px;
	friend class derivative_py;
	friend class derivative_pz;
	std::vector<celestial_body>bodies;
	double G;
	long long int n_iter;
	double time_step;
	double time = 0;
	long long int n_bodies;
	long long int curr_body;
	void start_soln() { //Reading for all n_bodies from file
		std::ifstream input_file("init_conds.txt");
		double temp;
		for (long long int o = 0; o < n_bodies; o++) {
			celestial_body new_body;
			input_file >> (new_body).mass;
			for (int m = 0; m < 3; m++) {
				input_file >> temp;
				(new_body.pos).push_back(temp);
				(new_body.provisional_pos).push_back(temp);
			}

			for (int m = 0; m < 3; m++) {
				input_file >> temp;
				(new_body.mom).push_back(temp);
				(new_body.provisional_mom).push_back(temp);
			}
			bodies.push_back(new_body);
		}
		input_file.close();
		
	}
public:
	n_body_data(double a, long long int t, double dt, long long int n) {
		G = a;
		n_iter = t;
		time_step = dt;
		n_bodies = n;
		start_soln();

	}
};
//Classes that will produce derivatives using operator overloading. Called by rk4 class.
class derivative_x {
	double y;
	double z;
	n_body_data& data;
	double x_dot = 0;

public:
	derivative_x(double y1, double z1, n_body_data& data1) : y(y1), z(z1), data(data1) {}

	double operator()(double x, double t) {
		
		for (long long int p = 0; p < data.n_bodies; p++) { //Force due to all n bodies taken into account
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt(((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2])));
			x_dot = x_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (x - data.bodies[p].pos[0]) *(denom * denom * denom);
			
		}

		return data.bodies[data.curr_body].mom[0] + x_dot * (t - data.time); //Euler method for finding momentum at
		// points within interval(in order to find derivatives of position for rk4). More  precise value of momentum only at
		// endpoints of each interval, using rk4 with the derivative_px, derivative_py, derivative_pz classes below.

	}
};

class derivative_y {
	double x;
	double z;
	n_body_data& data;
	double y_dot = 0;

public:
	derivative_y(double x1, double z1, n_body_data& data1) : x(x1), z(z1), data(data1) {}

	double operator()(double y, double t) {
		

		for (long long int p = 0; p < data.n_bodies; p++) {
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt(((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2])));
			y_dot = y_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (y - data.bodies[p].pos[1]) * (denom * denom * denom);


		}

		return data.bodies[data.curr_body].mom[1] + y_dot * (t - data.time);


	}
};

class derivative_z {
	double x;
	double y;
	n_body_data& data;
	double z_dot = 0;

public:
	derivative_z(double x1, double y1, n_body_data& data1) : x(x1), y(y1), data(data1) {}

	double operator()(double z, double t) {
		

		for (long long int p = 0; p < data.n_bodies; p++) {
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt(((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2])));
			z_dot = z_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (z - data.bodies[p].pos[2]) * (denom * denom * denom);


		}

		return data.bodies[data.curr_body].mom[2] + z_dot * (t - data.time);


	}
};

class derivative_px {
	double y;
	double z;
	n_body_data& data;
	double px_dot = 0;

public:
	derivative_px(double y1, double z1, n_body_data& data1) : y(y1), z(z1), data(data1) {}

	double operator()(double x, double t) {
		

		for (long long int p = 0; p < data.n_bodies; p++) {
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt(((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2])));
			px_dot = px_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (x - data.bodies[p].pos[0]) * (denom * denom * denom);


		}

		return px_dot;


	}
};

class derivative_py {
	double x;
	double z;
	n_body_data& data;
	double py_dot = 0;

public:
	derivative_py(double x1, double z1, n_body_data& data1) : x(x1), z(z1), data(data1) {}

	double operator()(double y, double t) {
		

		for (long long int p = 0; p < data.n_bodies; p++) {
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt(((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2])));
			py_dot = py_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (y - data.bodies[p].pos[1]) * (denom * denom * denom);


		}

		return py_dot;


	}
};

class derivative_pz {
	double x;
	double y;
	n_body_data& data;
	double pz_dot = 0;

public:
	derivative_pz(double x1, double y1, n_body_data& data1) : x(x1), y(y1), data(data1) {}

	double operator()(double z, double t) {
		

		for (long long int p = 0; p < data.n_bodies; p++) {
			if (p==data.curr_body) {
				continue;
			}
			double denom = find_fast_sqrt((x - data.bodies[p].pos[0]) * (x - data.bodies[p].pos[0]) + (y - data.bodies[p].pos[1]) * (y - data.bodies[p].pos[1]) + (z - data.bodies[p].pos[2]) * (z - data.bodies[p].pos[2]));
			
			pz_dot = pz_dot - data.G * data.bodies[p].mass * data.bodies[data.curr_body].mass * (z - data.bodies[p].pos[2]) * (denom * denom * denom);

		}

		return pz_dot;


	}
};



class implement_runge_kutta_step {
	double var;
	double t;
	double delta_t;
	double d1 = 0;
	double d2 = 0;
	double d3 = 0;
	double d4 = 0;
	double final_answer = 0;
public:
	implement_runge_kutta_step(double time_step, double var1, double t1) {
		delta_t = time_step;
		var = var1;
		t = t1;
		
	}

	template<class derivative> //Generic input giving derivatives
	double operator()(derivative d) {

		d1 = d(var, t);
		d2 = d(var + d1 * delta_t / 2, t + delta_t / 2);
		d3 = d(var + d2 * delta_t / 2, t + delta_t / 2);
		d4 = d(var + d3 * delta_t, t + delta_t);
		final_answer = (d1 + 2 * d2 + 2 * d3 + d4) / 6;
		return final_answer;
	}

};

class n_body_solution {
public:
	n_body_solution(n_body_data& data_input) {
		n_body_data& data1= data_input;
		for (long long int i = 0; i < data1.n_iter; i++) { //For each iteration, goes over every body, takes into account forces from all other bodies
			for (long long int j = 0; j < data1.n_bodies; j++) {

				data1.curr_body = j;

				implement_runge_kutta_step rx(data1.time_step, data1.bodies[j].pos[0], data1.time); //Create rk4 instance
				implement_runge_kutta_step ry(data1.time_step, data1.bodies[j].pos[1], data1.time);
				implement_runge_kutta_step rz(data1.time_step, data1.bodies[j].pos[2], data1.time);
				//To update values for all bodies in lockstep, updated values stored in provisional_pos, provisional_mom,
				//and pos, mom are updated after the loop runs through all bodies.
				data1.bodies[j].provisional_pos[0] = data1.bodies[j].provisional_pos[0] + rx(derivative_x(data1.bodies[j].pos[1], data1.bodies[j].pos[2], data1)) * data1.time_step / data1.bodies[j].mass; //Pass the derivative function for rk4 to be implemented
				data1.bodies[j].provisional_pos[1] = data1.bodies[j].provisional_pos[1] + ry(derivative_y(data1.bodies[j].pos[0], data1.bodies[j].pos[2], data1)) * data1.time_step / data1.bodies[j].mass;
				data1.bodies[j].provisional_pos[2] = data1.bodies[j].provisional_pos[2] + rz(derivative_z(data1.bodies[j].pos[0], data1.bodies[j].pos[1], data1)) * data1.time_step / data1.bodies[j].mass;
				
				data1.bodies[j].pos_history.push_back(data1.bodies[j].pos);

				implement_runge_kutta_step rpx(data1.time_step, data1.bodies[j].pos[0], data1.time);
				implement_runge_kutta_step rpy(data1.time_step, data1.bodies[j].pos[1], data1.time);
				implement_runge_kutta_step rpz(data1.time_step, data1.bodies[j].pos[2], data1.time);

				data1.bodies[j].provisional_mom[0] = data1.bodies[j].provisional_mom[0] + rpx(derivative_px(data1.bodies[j].pos[1], data1.bodies[j].pos[2], data1)) * data1.time_step;
				data1.bodies[j].provisional_mom[1] = data1.bodies[j].provisional_mom[1] + rpy(derivative_py(data1.bodies[j].pos[0], data1.bodies[j].pos[2], data1)) * data1.time_step;
				data1.bodies[j].provisional_mom[2] = data1.bodies[j].provisional_mom[2] + rpz(derivative_pz(data1.bodies[j].pos[0], data1.bodies[j].pos[1], data1)) * data1.time_step;

				data1.bodies[j].mom_history.push_back(data1.bodies[j].mom);

			}

			data1.time = data1.time + data1.time_step;  //Update time
		    for (long long int m = 0; m < data1.n_bodies; m++) {
				data1.bodies[m].pos = data1.bodies[m].provisional_pos;
				data1.bodies[m].mom = data1.bodies[m].provisional_mom;
		    }

		}
		//Send output to text files, with space between body entries
		std::ofstream output_file_1("data_pos.txt");
		for (long long int l = 0; l < data1.n_bodies; l++) {
			for (long long int k = 0; k < data1.n_iter; k++) {
				output_file_1 << data1.bodies[l].pos_history[k][0] << " " << data1.bodies[l].pos_history[k][1] << " " << data1.bodies[l].pos_history[k][2] << "\n";
			}
			output_file_1 << "\n";
		}
		output_file_1.close();

		std::ofstream output_file_2("data_mom.txt");
		for (long long int l = 0; l < data1.n_bodies; l++) {
			for (long long int k = 0; k < data1.n_iter; k++) {
				output_file_2 << data1.bodies[l].mom_history[k][0] << " " << data1.bodies[l].mom_history[k][1] << " " << data1.bodies[l].mom_history[k][2] << "\n";
			}
			output_file_2 << "\n";
		}
		output_file_2.close();


	}

};
int main(void) {
	n_body_data data1(1,50000,1E-16,2); //Relevant data for value of G, no. of iterations, time step, all held in data1.
	n_body_solution sol1(data1); //Implementation of solution

	return 0;
}