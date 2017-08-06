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
	
	//cout << "before  initialization for loop" << endl;
	
	for (int i=0; i < num_particles; i++){
		particles.push_back(Particle());
		particles[i].id = i;
		particles[i].weight = 1;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_psi(gen);
	}
	
	//cout << "after initialization for loop" << endl;
	//cout << "particles size" << particles.size() << endl;
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
	
	for (int i=0; i < num_particles; i++){
		particles[i].x += dist_v(gen)/dist_yr(gen)*(sin(particles[i].theta + dist_yr(gen)*delta_t)-sin(particles[i].theta));
		particles[i].y += dist_v(gen)/dist_yr(gen)*(cos(particles[i].theta) - cos(particles[i].theta + dist_yr(gen)*delta_t));
		particles[i].theta += dist_yr(gen)*delta_t;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//the first parameter has a rather misleading name: it is just the list of landmarks as provided by the map
	
	double min_distance = 10000;
	
	for (int i=0; i < num_particles; i++){ //loop through each particle, 
		for (int j=0; j < observations.size(); j++){ // for each particle, loop through each observations, and convert observation to map coordinates 
			
			// convert sensor observations to map coordinates
			//cout << "convert sensor observations to map coordinates" << endl;
			observations[j].x = particles[i].x*cos(particles[i].theta) - particles[i].y*sin(particles[i].theta) + observations[j].x;
			observations[j].y = particles[i].x*sin(particles[i].theta) - particles[i].y*cos(particles[i].theta) + observations[j].y;
			//cout << "convert sensor observations to map coordinates end" << endl;
			
			for (int k=0; k < predicted.size(); k++){ // loop through each landmark, and calculate distance to current observation in map coordinates
				
				//Calculate distance between landmark and observatin (in map coordinates)
				double distance = dist(predicted[k].x, predicted[k].y, observations[j].x, observations[j].y);  
								
				if (distance < min_distance) {
					min_distance = distance;
					
					//cout << "before assigining new values" << endl;
					//cout << "i: " << i << " type " << typeid(i).name() << endl;
					//cout << "obsx " << typeid(observations[j].x ).name() << endl;
					//cout << "obsy " << observations[j].x << endl;
					
					//cout << "left side" << particles[i].associations[0] << endl;
					
					//particles[i].associations[0] = 1; //predicted[j].id;
					//cout << "past1" << endl;
					vector<double> obsxx = {observations[j].x};
					vector<double> obsyy = {observations[j].y};
					particles[i].sense_x = obsxx; //landmark
					//cout << "past2" << endl;
					particles[i].sense_y = obsyy;
					//cout << "past3" << endl;
					particles[i].id = predicted[k].id; 
				}
			}
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
	
	//convert ‘std::vector<Map::single_landmark_s>’ to ‘std::vector<LandmarkObs>’
	
	std::vector<LandmarkObs> landmarkData;
	//cout << "data type conversion begin" << endl;
	for (int i=0; i < map_landmarks.landmark_list.size(); i++){
			landmarkData.push_back(LandmarkObs());
			landmarkData[i].id = map_landmarks.landmark_list[i].id_i;
			landmarkData[i].x = map_landmarks.landmark_list[i].x_f;
			landmarkData[i].y = map_landmarks.landmark_list[i].y_f;
	}
	//cout << "data type conversion end, data association begin" << endl;
	
	ParticleFilter::dataAssociation(landmarkData, observations);
	
	for (int i=0; i < num_particles; i++) { //loop through all the particles
		
		for (int j=0; j < particles[i].associations.size(); j++) { //each particle is associated with a number of observations, so loop through all observations
		// these associated observations are converted to map coordinates (e.g. GPS)
			particles[i].weight = 1 / (2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-(pow((particles[i].sense_x[j] - map_landmarks.landmark_list[particles[i].id].x_f ),2)/(2*pow(std_landmark[0],2)) + pow(( particles[i].sense_y[j] - map_landmarks.landmark_list[particles[i].id].y_f ),2 )/(2*pow(std_landmark[1],2))   ));
		}
		
		
	}
        
	
	
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	int N = num_particles;
	vector<Particle> new_particles;
	int index = rand() % N;
	double beta = 0.0;
	
	// find largest weight
	decltype(particles)::iterator ParticleIterator;
	
	ParticleIterator = max_element(begin(particles), end(particles),
		[] (Particle const& p1, Particle  const& p2)
	    {
	        return p1.weight < p2.weight;
	    });
	
	double mw = particles[std::distance(begin(particles), ParticleIterator)].weight;
	//cout << "here!" << endl;
	//double mw = particles[3].weight ;
	for (int i=0; i < N; i++){
		//cout << "1" << endl;
		beta = beta + rand() % 3 * mw;
		//cout << "2" << endl;
		while (particles[index].weight < beta) {
			//cout << "3" << endl;
			beta = beta - particles[index].weight;
			//cout << "4" << endl;
			index = index + 1;
			//cout << "5" << endl;
		new_particles.push_back(Particle());
		new_particles[i] = particles[index];
		//cout << "6" << endl;
		}
	}
	
	
	particles = new_particles;
	//cout << "7" << endl;
	
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
