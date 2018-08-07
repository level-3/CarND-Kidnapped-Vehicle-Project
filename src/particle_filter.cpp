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

	return;
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
			if (fabs(yaw_rate) < 0.00001)
			{
			particles[i].x += velocity * delta_t * cos(theta) ;
			particles[i].y += velocity * delta_t * sin(theta) ;			
			}
			else
			{
			particles[i].x += (velocity / yaw_rate) * ( sin(theta + yaw_rate * delta_t) - sin ( theta ) ) ;
			particles[i].y += (velocity / yaw_rate) * ( cos(theta) - cos ( theta + yaw_rate * delta_t ) ) ;
			particles[i].theta += theta + yaw_rate * delta_t;
			}		


			//add noise
			default_random_engine gen;
			particles[i].x = dist_x(gen);
			particles[i].y = dist_y(gen);
			particles[i].theta = dist_theta(gen);

		};

	return;




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


	return;
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

	for (int i = 0; i < num_particles ; i++)
	{
		double part_x = particles[i].x ;
		double part_y = particles[i].y ; 
		double part_theta = particles[i].theta ;

		vector<LandmarkObs> predictions;

		std::vector<Map::single_landmark_s> land_list = map_landmarks.landmark_list;

		for ( unsigned int j = 0; j < land_list.size() ; j++)
		{
			int land_id = land_list[j].id_i;			
			float land_x = land_list[j].x_f;
			float land_y = land_list[j].y_f;

			float dX = part_x - land_x;
			float dY = part_y - land_y;

			if (dX*dX + dY*dY <= sensor_range*sensor_range)
			{
				LandmarkObs pred = {land_id,land_x,land_y};
				predictions.push_back( pred );	
			}
			
		}

		// transform obs coord
		vector<LandmarkObs> transformed_obs;

		for ( unsigned int j = 0 ; j < observations.size() ; j++)
		{
			double t_x = cos(part_theta) * observations[j].x - sin(part_theta) * observations[j].y + part_x ;
			double t_y = cos(part_theta) * observations[j].x - cos(part_theta) * observations[j].y + part_y ;

			transformed_obs.push_back( LandmarkObs{observations[j].id,t_x,t_y} );
		}

		dataAssociation(predictions, transformed_obs);
 
		// reset weight
		particles[i].weight = 1.0;

		for ( unsigned int j = 0 ; j < transformed_obs.size() ; j++)
		{
			double pred_x, pred_y;			
			double obs_x = transformed_obs[j].x;
			double obs_y = transformed_obs[j].y;	
			int id_landmark = transformed_obs[j].id;

			for (unsigned int k = 0 ; k <predictions.size();k++)
			{
				if ( id_landmark == predictions[k].id )
				{
					pred_x = predictions[k].x;
					pred_y = predictions[k].y;					
				}
			}

		// calculate weights with MV Gaussian

		double s_x = std_landmark[0];
		double s_y = std_landmark[1];

		double gauss_norm = (1/(2 * M_PI * s_x*s_y ) );
		double exponent = - pow( pred_x - obs_x , 2 ) / ( 2 * pow(s_x, 2) ) + pow(pred_y-obs_y,2) / ( 2* pow(s_y, 2) );
		double weight = gauss_norm * exponent;
		
		particles[i].weight *= weight;

		}

	}

	return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
return;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
