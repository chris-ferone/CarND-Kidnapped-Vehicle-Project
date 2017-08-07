/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <typeinfo>
#include <iomanip> 

#include "particle_filter.h"
#include <cmath> 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
		
	// This line create normal (Gaussian) distributions for x, y, and psi
	
	num_particles = 10;
	
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);
	default_random_engine gen;
	
	cout << "gps x" << x << endl;
	cout << "gps y" << y << endl;
	cout << "gps theta" << theta << endl;
	
	for (int i=0; i < num_particles; i++){
		particles.push_back(Particle());
		particles[i].id = i;
		particles[i].weight = 1.0;
		particles[i].x = dist_x(gen);
		cout << "x initial: " << particles[i].x << endl;
		particles[i].y = dist_y(gen);
		cout << "y initial: " << particles[i].y << endl;
		particles[i].theta = dist_psi(gen);
		cout << "theta initial: " << particles[i].theta << endl;
	}
	
	is_initialized = true;
	


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	
	default_random_engine gen;
	
	for (int i=0; i < num_particles; i++){
		// Calculate the mean
		if(abs(yaw_rate)<0.0001) {
 			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta); }
		else {
			particles[i].x += velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta)); 
			particles[i].y += velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate*delta_t;}
		//create a normal distribution to represent 0 mean noise with standard deviation sigma
		normal_distribution<double> dist_x(0 , std_pos[0]);
		normal_distribution<double> dist_y(0 , std_pos[1]);
		normal_distribution<double> dist_psi(0 , std_pos[2]);
		//sample from normal distribution and add to mean
		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_psi(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	int num_observations = observations.size();
	vector<int> assoc(num_observations); //Each particle has a some number of landmarks associated with (depends on what was in sensor's field of view)
	vector<double> obsxx(num_observations);
	vector<double> obsyy(num_observations);
		
	
	for (int i=0; i < num_particles; i++){ //loop through each particle, 
		for (int j=0; j < observations.size(); j++){ // for each particle, loop through each observations, and convert observation to map coordinates 
			// convert sensor observations to map coordinates

			LandmarkObs observMAP;
			observMAP.x = 0;
			observMAP.y = 0;

			observMAP.x = observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			observMAP.y = observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y;

			obsxx[j] = observMAP.x;
			obsyy[j] = observMAP.y;

			int min_distance_k;
			double min_distance = 10000;
			for (int k=0; k < map_landmarks.landmark_list.size(); k++){ // loop through each landmark, and calculate distance to current observation in map coordinates
				
				//Calculate distance between landmark and observatin (in map coordinates)
				double distance = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, observMAP.x, observMAP.y);  		
				if (distance < min_distance) {
					min_distance = distance;
					min_distance_k = k;

				}
			}
			assoc[j] = min_distance_k+1;
		}
		particles[i].sense_x = obsxx;
		particles[i].sense_y = obsyy;
		particles[i].associations = assoc;
	}
	
	weights.clear();
	
	for (int i=0; i < num_particles; i++) { //loop through all the particles
		particles[i].weight = 1.0;
			for (int j=0; j < observations.size(); j++){  
				particles[i].weight *= 1.0 / (2.0*M_PI*std_landmark[0]*std_landmark[1]) * exp(-( 
					pow( (particles[i].sense_x[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].x_f ) ,2) / (2.0*pow(std_landmark[0],2)) + 
					pow( (particles[i].sense_y[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].y_f ) ,2) / (2.0*pow(std_landmark[1],2)) ));
			}
			
			weights.push_back(particles[i].weight);
		}
			
		
	
        
	
	
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> new_particles(num_particles);
	default_random_engine generator;
	
	discrete_distribution<> dist (weights.begin(), weights.end());
	
	for (int i = 0; i < num_particles; i++)
	{
		int idx = dist(generator);
		new_particles[i] = particles[idx];
	}

	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
