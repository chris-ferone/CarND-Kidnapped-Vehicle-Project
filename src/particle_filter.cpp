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
	
	num_particles = 2;
	
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);
	default_random_engine gen;
	
	cout << "gps x" << x << endl;
	cout << "gps y" << y << endl;
	cout << "gps theta" << theta << endl;
	
	//cout << "before  initialization for loop" << endl;
	
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
		cout << "before prediciton: " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << " " << endl;
		particles[i].x += dist_v(gen)/dist_yr(gen)*(sin(particles[i].theta + dist_yr(gen)*delta_t)-sin(particles[i].theta));
		particles[i].y += dist_v(gen)/dist_yr(gen)*(cos(particles[i].theta) - cos(particles[i].theta + dist_yr(gen)*delta_t));
		particles[i].theta += dist_yr(gen)*delta_t;
		cout << "after prediciton: " << particles[i].x << " " << particles[i].y << " " << particles[i].theta << " " << endl;
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
	
	//convert ‘std::vector<Map::single_landmark_s>’ to ‘std::vector<LandmarkObs>’
	
/* 	std::vector<LandmarkObs> landmarkData;
	cout << "data type conversion begin" << endl;
	for (int i=0; i < map_landmarks.landmark_list.size(); i++){
			landmarkData.push_back(LandmarkObs());
			landmarkData[i].id = map_landmarks.landmark_list[i].id_i;
			landmarkData[i].x = map_landmarks.landmark_list[i].x_f;
			landmarkData[i].y = map_landmarks.landmark_list[i].y_f;
	}
	cout << "data type conversion end, data association begin" << endl;
	
	ParticleFilter::dataAssociation(landmarkData, observations); */
	
/* 	cout << "map index 0 x" << map_landmarks.landmark_list[0].x_f << endl;
	cout << "map index 0 x" << map_landmarks.landmark_list[0].y_f << endl;
	cout << "map index 0 x" << map_landmarks.landmark_list[0].id_i << endl;
	
	cout << "map index 1 x" << map_landmarks.landmark_list[1].x_f << endl;
	cout << "map index 1 x" << map_landmarks.landmark_list[1].y_f << endl;
	cout << "map index 1 x" << map_landmarks.landmark_list[1].id_i << endl; */
	
	int num_observations = observations.size();
	vector<int> assoc(num_observations); //Each particle has a some number of landmarks associated with (depends on what was in sensor's field of view)
	vector<double> obsxx(num_observations);
	vector<double> obsyy(num_observations);
		
	
	for (int i=0; i < num_particles; i++){ //loop through each particle, 
	

		//cout << "i: " << i <<endl;
		for (int j=0; j < observations.size(); j++){ // for each particle, loop through each observations, and convert observation to map coordinates 
			
			//cout << "observation x: " << observations[j].x << " || observation y: " << observations[j].y << endl;
			//cout << "particle x: " << particles[i].x << " || particle y: " << particles[i].y <<  " || particle theta: "   << particles[i].theta << endl;
			
			// convert sensor observations to map coordinates
			//cout << "convert sensor observations to map coordinates" << endl;
			LandmarkObs observMAP;
			observMAP.x = 0;
			observMAP.y = 0;
			//observMAP.x = particles[i].x*cos(particles[i].theta) - particles[i].y*sin(particles[i].theta) + observations[j].x;
			//observMAP.y = particles[i].x*sin(particles[i].theta) + particles[i].y*cos(particles[i].theta) + observations[j].y;
			observMAP.x = observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			observMAP.y = observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y;
			
			//cout << "convert sensor observations to map coordinates end" << endl;
			
			obsxx[j] = observMAP.x;
			obsyy[j] = observMAP.y;
			
			//cout << "observation x MAP: " << observMAP.x << " || observation y MAP: " << observMAP.y << endl;
			
			//cout << "observation x MAP: " << obsxx[j] << " || observation y MAP: " << obsyy[j] << endl;
			
			int min_distance_k;
			double min_distance = 10000;
			for (int k=0; k < map_landmarks.landmark_list.size(); k++){ // loop through each landmark, and calculate distance to current observation in map coordinates
				
				//Calculate distance between landmark and observatin (in map coordinates)
				double distance = dist(map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, observMAP.x, observMAP.y);  
				//cout << "distance " << distance <<endl;			
				if (distance < min_distance) {
					min_distance = distance;
					min_distance_k = k;

				}
			}
			assoc[j] = min_distance_k+1;
			//cout << "minimum distance " << min_distance <<endl;	
			//cout << "associated landmark ID " << min_distance_k + 1<<endl;
		}
		particles[i].sense_x = obsxx;
		particles[i].sense_y = obsyy;
		particles[i].associations = assoc;
	}
	
	
	
	for (int i=0; i < num_particles; i++) { //loop through all the particles
		particles[i].weight = 1.0;
		//for (int j=0; j < particles[i].associations.size(); j++) { //each particle is associated with a number of observations, so loop through all observations
		// these associated observations are converted to map coordinates (e.g. GPS)
/* 			cout << "sense_x: " << fixed << setprecision(12) << particles[i].sense_x[0] << endl;
			cout << "land_x: " << fixed << setprecision(12) << map_landmarks.landmark_list[particles[i].id+1].x_f << endl;
			cout << "sense_y: " << fixed << setprecision(12) << particles[i].sense_y[0] << endl;
			cout << "land_y: " << fixed << setprecision(12) << map_landmarks.landmark_list[particles[i].id+1].y_f << endl; */
			
				
			//double bb = pow((particles[i].sense_x[j] - map_landmarks.landmark_list[particles[i].id+1].x_f ),2)/(2*pow(std_landmark[0],2));
			//cout << "bb: " << fixed << setprecision(12) << bb << endl;
			//double cc = pow((particles[i].sense_y[j] - map_landmarks.landmark_list[particles[i].id+1].y_f ),2)/(2*pow(std_landmark[1],2));
			//cout << "cc: " << fixed << setprecision(12) << cc << endl;
			//double aa = exp(-( bb + cc));
			//cout << "aa: " << fixed << setprecision(12) << aa << endl;
			for (int j=0; j < observations.size(); j++){
 				 //  cout << "sensed x: " << particles[i].sense_x[j] << " actual x: " << map_landmarks.landmark_list[particles[i].associations[j]-1].x_f << endl;
				//cout << "sensed y: " << particles[i].sense_y[j] << " actual y: " << map_landmarks.landmark_list[particles[i].associations[j]-1].y_f << endl;
				//cout << "x diff: " << (particles[i].sense_x[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].x_f ) << 
				//"y diff: " << (particles[i].sense_y[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].y_f ) << endl;    
				particles[i].weight *= 1.0 / (2.0*M_PI*std_landmark[0]*std_landmark[1]) * exp(-( 
					pow( (particles[i].sense_x[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].x_f ) ,2) / (2.0*pow(std_landmark[0],2)) + 
					pow( (particles[i].sense_y[j] - map_landmarks.landmark_list[particles[i].associations[j]-1].y_f ) ,2) / (2.0*pow(std_landmark[1],2)) ));
				//cout << "pre particle weights: " << fixed << setprecision(20) << particles[i].weight << endl;
			}
			
			//cout << "updated particle weights: " << fixed << setprecision(20) << particles[i].weight << endl;
			weights.push_back(particles[i].weight);
		}
			
		
	
        
	
	
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> new_particles(num_particles);
/* 	int index = rand() % num_particles;
	double beta = 0.0;
	
	// find largest weight
	decltype(particles)::iterator ParticleIterator;
	
	ParticleIterator = max_element(begin(particles), end(particles),
		[] (Particle const& p1, Particle  const& p2)
	    {
	        return p1.weight < p2.weight;
	    });
	
	double mw = particles[std::distance(begin(particles), ParticleIterator)].weight;
	cout << "mw " << mw << endl;
	//double mw = particles[3].weight ;
	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	for (int i=0; i < num_particles; i++){
		//cout << "1" << endl;

		uniform_real_distribution<> dis(0, 2*mw);
		beta = beta + dis(gen);
		cout << "beta2: " << beta << endl;
		//cout << "weight: " << particles[index].weight << endl;
		while (particles[index].weight < beta) {
			cout << "beta3: " << beta << endl;
			beta = beta - particles[index].weight;
			cout << "beta4: " << beta << endl;
			index = index + 1;
			cout << "5" << endl;
		}
		//new_particles.push_back(Particle());
		new_particles[i] = particles[index];
		cout << "new particle added: " << i << endl;
		
	} */
	//cout << "particles size" << particles.size() << " num_particles " << num_particles << endl;
	default_random_engine generator;
	
	discrete_distribution<> distribution (0, num_particles-1);

	for (int i = 0; i < num_particles; i++)
	{
		int idx = distribution(generator);
		new_particles[i] = particles[idx];
		//cout <<"i: " << i << endl;
	}
	
	
	particles = new_particles;
	//cout << "7" << endl;
	for (int i = 0; i < particles.size(); i++)
	{
		cout <<"i: " << i << "   x final: " << particles[i].x << endl;
	}
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
