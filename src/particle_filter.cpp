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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	// STUDENT'S NOTE: Should the added random Gaussian noise be different from the uncertainty
	// already accounted for through GPS std? This is ambiguous. See also my note in prediction().
	
	// Prepare normal distributions
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// Create particles and initialize them with (soft) GPS position
	for(int ii = 0; ii < this->num_particles; ii = ii + 1) {
		Particle p = Particle();
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1.0;

		this->particles.push_back(p);		
	}
	// Signal that the filter has been initialized
	this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.

	// STUDENT'S NOTE: Again, we do not add measurements and add random Gaussian noise. What we do is to
	// draw the measurements from a normal distribution with a certain mean (the measurement itself) and (co)variance.
	// Even though programmatically those two concepts are the same they shouldn't (in my opinion) be explained like that.

	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Distributions for motion uncertainty
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	// Our effective "zero" value
	double EPSILON = 1e-6;

	// Update (x, y, theta) for each particle (checking if yaw_rate == 0)
	for(int ii = 0; ii < this->num_particles; ii = ii + 1) {
		double theta = this->particles[ii].theta;
		double update;

		if(abs(yaw_rate) > EPSILON) {
			update = velocity/yaw_rate * (sin(theta + yaw_rate*delta_t) - sin(theta));
			this->particles[ii].x = this->particles[ii].x + update + dist_x(gen);

			update = velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate*delta_t));
			this->particles[ii].y = this->particles[ii].y + update + dist_y(gen);

			this->particles[ii].theta = this->particles[ii].theta + yaw_rate*delta_t + dist_theta(gen);
		} else {
			update = velocity*delta_t*cos(theta);
			this->particles[ii].x = this->particles[ii].x + update;

			update = velocity*delta_t*sin(theta);
			this->particles[ii].y = this->particles[ii].y + update;
			
		}
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.

	// STUDENT NOTE: The TODO doesn't make any sense. E.g. who is "this particular landmark"? The
	// header file is not helping either.

	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::findNearestNeighbours(std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Perform data association by nearest neighbour approach.

	// For each observation
	for(unsigned int ii = 0; ii < observations.size(); ii = ii + 1) {		
		double ref_x = observations[ii].x;
		double ref_y = observations[ii].y;
		
		// ...find the closest match (nearest landmark)
		// By definition each observation has one.
		double nn_distance = 1e10;
		int nn_index = 0;
		for(unsigned int jj = 0; jj < map_landmarks.landmark_list.size(); jj = jj + 1) {
			double x = map_landmarks.landmark_list[jj].x_f;
			double y = map_landmarks.landmark_list[jj].y_f;
			double index = map_landmarks.landmark_list[jj].id_i;

			// Compute squared delta in each axis and the resulting Euclidean distance
			double delta_x2 = (ref_x-x)*(ref_x-x);
			double delta_y2 = (ref_y-y)*(ref_y-y);
			double distance = sqrt(delta_x2 + delta_y2);

			// Update if distance is smaller than the current best
			if(distance < nn_distance) {
				nn_distance = distance;
				nn_index = index;
			}
		}

		// Cross check
		// In order to avoid duplicities, we check that, given the observation,
		// the nearest neighbour for its closest landmark is the observation itself.
		bool cross_check_ok = true;
		// Given the selected landmark
		ref_x = map_landmarks.landmark_list[nn_index].x_f;
		ref_y = map_landmarks.landmark_list[nn_index].y_f;
		// Look for its closest observation
		for(unsigned int jj = 0; jj < observations.size(); jj = jj + 1) {
			double x = observations[ii].x;
			double y = observations[ii].y;
			double delta_x2 = (ref_x-x)*(ref_x-x);
			double delta_y2 = (ref_y-y)*(ref_y-y);
			double distance = sqrt(delta_x2 + delta_y2);			
			if(distance < nn_distance) {
				// Cross-check failed, the nearest landmark for this observation and
				// the nearest observation for this landmark are not the same.
				cross_check_ok = false;
				break;
			}
		}

		if(cross_check_ok) {
			observations[ii].id = nn_index;
		} else {
			observations[ii].id = 0;
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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

	// Transform the observation as seen from each particle and
	// express it in the global (map) coordinate system
	for(int ii = 0; ii < this->num_particles; ii = ii + 1) {
		double xp = this->particles[ii].x;
		double yp = this->particles[ii].y;
		double theta = this->particles[ii].theta;

		std::vector<LandmarkObs> predicted;
		for(unsigned int jj = 0; jj < observations.size(); jj = jj + 1) {
			double xc = observations[jj].x;
			double yc = observations[jj].y;

			double distance = sqrt(xc*xc + yc*yc);
			if(distance > sensor_range) {
				cout << "Observation out of range: " << distance << endl;
				continue;
			}

			double xm = xp + (xc*cos(theta)) - (yc*sin(theta));
			double ym = yp + (xc*sin(theta)) + (yc*cos(theta));

			LandmarkObs prediction;
			prediction.x = xm;
			prediction.y = ym;
			predicted.push_back(prediction);
		}

		// Perform data association
		findNearestNeighbours(predicted, map_landmarks);

		// Compute weights
		double new_weight = 1.0;
		for(unsigned int kk = 0; kk < predicted.size(); kk = kk + 1) {
			// We have to substract 1 since id starts at 1 and not 0
			int match = predicted[kk].id - 1;
			if(match == 0) {
				continue;
			}
			double dx = predicted[kk].x - map_landmarks.landmark_list[match].x_f;
			double dy = predicted[kk].y - map_landmarks.landmark_list[match].y_f;

			if (sqrt(dx*dx + dy*dy) >= sensor_range) {
				continue;
			}			

			double expo = (dx*dx)/(2*std_landmark[0]*std_landmark[0]) + (dy*dy)/(2*std_landmark[1]*std_landmark[1]);
			new_weight = new_weight * exp(-expo) / sqrt(2*std_landmark[0]*std_landmark[1]*M_PI);			
		}
		this->particles[ii].weight = new_weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	// NOTE> We use std::discrete_distribution instead of resampling wheel.

	// Collect the weights
	std::vector<double> weights (this->num_particles);
	for(int ii = 0; ii < this->num_particles; ii = ii + 1) {
		weights[ii] = this->particles[ii].weight;
	}
	
	// Set up a discrete distribution according to the weights
	default_random_engine gen;
	std::discrete_distribution<> distribution(weights.begin(), weights.end());

	// Resampling process
	std::vector<Particle> resampled_particles;
	for(int ii = 0; ii < this->num_particles; ii = ii + 1) {
		int index = distribution(gen);
		resampled_particles.push_back(this->particles[index]);
	}
	this->particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
	//STUDENT'S NOTE: This function should return Particle object.
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
	// STUDENT'S NOTE: The variable "associations", as Particle attribute should be explained in Particle class
	// and not here. Idem for sense_x, sense_y.
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
	// STUDENT'S NOTE: In updateWeights() the instructions mention map's and vehicle's coordinate system.
	// Are world coordinates something else? I beleive they in fact refer to the map coordinates but why
	// are there two names for the same concept?

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;	
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
