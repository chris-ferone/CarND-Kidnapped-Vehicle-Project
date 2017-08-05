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

#include "particle_filter.h"
#include <cmath> 

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
		
	// This line create normal (Gaussian) distributions for x, y, and psi
	
	num_particles = 20;
	
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);
	default_random_engine gen;
	
	for (int i=0; i < size(num_particles); i++){
		Particle[i].id = i;
		Particle[i].Weight = 1;
		Particle[i].x = dist_x(gen);
		Particle[i].y = dist_y(gen);
		Particle[i].theta = dist_psi(gen);
	}
	
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	normal_distribution<double> dist_v(velocity, std_pos[0]);
	normal_distribution<double> dist_yr(yaw_rate, std_pos[1]);
	default_random_engine gen;
	
	for (int i=0; i < size(num_particles); i++){
		Particle[i].x += dist_v(gen)/dist_yr(gen)*(sin(Particle.theta + dist_yr(gen)*delta_t)-sin(Particle.theta));
		Particle[i].y += dist_v(gen)/dist_yr(gen)*(cos(Particle.theta) - cos(Particle.theta + dist_yr(gen)*delta_t));
		Particle[i].theta += dist_yr(gen)*detlta_t;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//the first parameter has a rather misleading name: it is just the list of landmarks as provided by the map
	
	
	// there are predicted measurements, observed measurements, and landmarks what do I assign to what?
	// dataAssociation is the process of associating sensor measurements with map landmarks
	double min_distance = 10000;
	
	for (int i=0; i < size(predicted); i++){ //loop through each observation, 
		for (int j=0; j < size(predicted); j++){ // for each observation, find the closest landmark; so loop through all landmarks, 
			
			// convert sensor observations to map coordinates
			observations.x = Particle[i].x*cos(Particle[i].theta) - Particle[i].y*sin(Particle[i].theta) + observation.x;
			observations.y = Particle[i].x*sin(Particle[i].theta) - Particle[i].y*cos(Particle[i].theta) + observation.y;
			
			distance = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);  
			if (distance < min_distance) {
				min_distance = distance;
				Particle[i].associations = Predicted[j].id_i;
				Particle[i].sense_x = observations.x; //landmark
				Particle[i].sense_y = observations.y;
			}
		
	}
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
	

	
	
	ParticleFilter::dataAssociation(map_landmarks, observations);
	
	for (int i=0; i < size(num_particles); i++) { //loop through all the particles
		
		for (int j=0; j < size(Particle[i].associations); j++) { //each particle is associated with a number of observations, so loop through all observations
		// these associated observations are converted to map coordinates (e.g. GPS)
			Particle[i].Weight = 1/(2*M_PI*sens_covar[0]*sens_covar[1])*exp(-(pow((Particle[i].sense_x - map_landmarks[Particle[i].id_i].x_f ),2)/(2*pow(sens_covar[0],2)) + pow(( Particle[i].sense_y - map_landmarks[Particle[i].id_i].y_f ),2 )/(2*pow(sens_covar[1],2))   ));
		}
		
		
	}
        
	
	
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	N = size(num_particles);
	std::vector<Particle> new_particles;
	int index = rand() % N;
	double beta = 0.0
	mw = Particles[max_element(std::begin(Particles), std::end(Particles))].weight ;
	for (int i=0; i < N; i++){
		beta = beta + rand() % 3 * mw;
		while Particle[index].weight < beta:
			beta = beta - Particle[index].weight;
			index = index + 1;
		new_particles[i] = Particle[index];
	}
	
	
	Particle = new_particles;
	
	
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
