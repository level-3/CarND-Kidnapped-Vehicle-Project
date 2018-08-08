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

	if (is_initialized == false) 
	{
		
		double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

		// Set standard deviations for x, y, and theta
		 std_x = std[0];
		 std_y = std[1];
		 std_theta = std[2];

		// create a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(x, std_x);

		// Create normal distributions for y and theta
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta, std_theta);

		// initialise number of particles
		num_particles  = 100;

		default_random_engine gen;

		for (int i=0; i<num_particles;i++)
		{
			Particle particle;
			particle.id = i;
			particle.x = dist_x(gen);
			particle.y = dist_y(gen);
			particle.theta = dist_theta(gen);	 
			particle.weight = 1.0;

			particles.push_back(particle);
		}

		is_initialized=true;
	}

}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
		
		double std_x, std_y, std_theta; 
		// get standard deviations for x, y, and theta
		 std_x = std_pos[0];
		 std_y = std_pos[1];
		 std_theta = std_pos[2];

		// create a normal (Gaussian) distribution for x
		normal_distribution<double> dist_x(0, std_x);
		normal_distribution<double> dist_y(0, std_y);
		normal_distribution<double> dist_theta(0, std_theta);

		for (int i = 0 ; i<num_particles;i++)
		{
			double theta = particles[i].theta;			
			if (fabs(yaw_rate) < 0.0001)
			{
				particles[i].x += velocity * delta_t * cos(theta) ;
				particles[i].y += velocity * delta_t * sin(theta) ;			
			}
			else
			{
				particles[i].x += velocity / yaw_rate * ( sin(theta + yaw_rate * delta_t) - sin ( theta ) ) ;
				particles[i].y += velocity / yaw_rate * ( cos(theta) - cos ( theta + yaw_rate * delta_t ) ) ;
				particles[i].theta +=  yaw_rate * delta_t;
			}		

			//add noise

			particles[i].x += dist_x(gen);
			particles[i].y += dist_y(gen);
			particles[i].theta += dist_theta(gen);

		};
}



void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < observations.size(); i++ )
	{
		LandmarkObs obs = observations[i];
		double min_dist = numeric_limits<double>::max();

		int map_id = -1;
		
		for (unsigned int j = 0; j < predicted.size(); j++)
		{
			LandmarkObs pred = predicted[j];

			double cur_dist = dist(obs.x,obs.y,pred.x,pred.y);
			
			if (cur_dist < min_dist)
			{
				min_dist = cur_dist ;
				map_id = pred.id;
			}

		}

		observations[i].id = map_id;
	}

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks)  {
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


  for (int i = 0; i < num_particles; i++) 
	{
		double part_x = particles[i].x;
		double part_y = particles[i].y;
		double theta = particles[i].theta;

    	// Find landmark range
		//double sensor_range_sqd = sensor_range * sensor_range;

		vector<LandmarkObs> land_inrange;

		for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) 
		{
			float land_x = map_landmarks.landmark_list[j].x_f;
			float land_y = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			
			double dX = part_x - land_x;
			double dY = part_y - land_y;

			if ( pow(dX,2) + pow(dY,2) <= pow(sensor_range,2) ) 
			{
				land_inrange.push_back(LandmarkObs{ id, land_x, land_y });
			}
		}

		// Transform observation coord
		vector<LandmarkObs> map_obs;

		for(unsigned int j = 0; j < observations.size(); j++) 
		{
			double x_m = cos(theta) * observations[j].x - sin(theta) * observations[j].y + part_x;
			double y_m = sin(theta) * observations[j].x + cos(theta) * observations[j].y + part_y;
			map_obs.push_back( LandmarkObs{ observations[j].id, x_m, y_m } );
		}

		dataAssociation(land_inrange, map_obs);

		particles[i].weight = 1.0;

		// Calc weights
		for(unsigned int j = 0; j < map_obs.size(); j++) 
			{
				double obs_x = map_obs[j].x;
				double obs_y = map_obs[j].y;
				int landmarkId = map_obs[j].id;
				double land_x, land_y;
				
				unsigned int k = 0;

				bool found = false;
				
				while( !found && k < land_inrange.size() ) 
				{
					if ( land_inrange[k].id == landmarkId)
					{
						found = true;
						land_x = land_inrange[k].x;
						land_y = land_inrange[k].y;
					}
					
				k++;
			}
			
			double std_range = std_landmark[0];
  			double std_bearing = std_landmark[1];

			double dX = obs_x - land_x;
			double dY = obs_y - land_y;
			
			double gauss_norm =  1 / ( 2 * M_PI * std_range * std_bearing );
			double exponent = exp( -( dX*dX / ( 2 * std_range*std_range ) + ( dY*dY / ( 2 * std_bearing*std_bearing ) ) ) );
			double weight = gauss_norm * exponent;
			
			if (weight == 0) 
			{
				particles[i].weight *= 0.0001;
			} 
			else 
			{
				particles[i].weight *= weight;
			}
		}
  	}
}



void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	

	vector<double> weights;  

	// set max weight
	double max_weight = numeric_limits<double>::min();

	for(int i = 0; i < num_particles; i++) 
	{
		weights.push_back(particles[i].weight);

		if ( particles[i].weight > max_weight ) 
		{
			max_weight = particles[i].weight;
		}
  	}

	// Creating distribution
	uniform_real_distribution<double> dist_double(0.0, max_weight);
	uniform_int_distribution<int> dist_int(0, num_particles - 1);

	// Create index
	int index = dist_int(gen);

	double beta = 0.0;

	// wheel
	vector<Particle> resampled_part;

	for(int i = 0; i < num_particles; i++) 
	{
		beta += dist_double(gen) * 2.0;

		while( beta > weights[index] ) 
		{
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		resampled_part.push_back(particles[index]);
	}

	particles = resampled_part;
}



Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
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
