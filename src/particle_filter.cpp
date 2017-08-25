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

    default_random_engine gen;

    //This number is determined by experiments
    num_particles = 200;

    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i =0; i < num_particles; i++)
    {

        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;

        particles.push_back(particle);

    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate){
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    normal_distribution<double> dist_x(0, std_pos[0]);
    normal_distribution<double> dist_y(0, std_pos[1]);
    normal_distribution<double> dist_theta(0, std_pos[2]);

    double x0,y0,theta0;

    //Check if the yaw rate is not a small value to avoid dividing by zero
    if (fabs(yaw_rate) > 0.00001)
    {
        for (int i =0; i < num_particles; i++)
        {
            x0 = particles[i].x;
            y0 = particles[i].y;
            theta0 = particles[i].theta;

            particles[i].x = x0 + (velocity/yaw_rate) * (sin(theta0 + yaw_rate*delta_t) - sin(theta0)) + dist_x(gen);
            particles[i].y = y0 + (velocity/yaw_rate) * (cos(theta0) - cos(theta0 + yaw_rate*delta_t)) + dist_y(gen);
            particles[i].theta = theta0 + yaw_rate*delta_t + dist_theta(gen);

        }

    }
    else
    {
        for (int i =0; i < num_particles; i++)
        {
            x0 = particles[i].x;
            y0 = particles[i].y;
            theta0 = particles[i].theta;

            particles[i].x = x0 + velocity * delta_t * cos(theta0) + dist_x(gen);
            particles[i].y = y0 + velocity * delta_t * sin(theta0) + dist_y(gen);
            particles[i].theta = theta0 + dist_theta(gen);

        }

    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    //loop over all the observed measurements
    for(int i =0; i<observations.size(); i++)
    {
        LandmarkObs curr_obs = observations[i];

        // initialize the minimum distance with the biggest double number
        double min_dist = std::numeric_limits<double>::max();

        //loop over all the predicted measurements
        for(int j =0; j < predicted.size(); j++)
        {

            //Calculate the distance between the predicted measurement and the current observation
            double dist_to_obs = dist(predicted[j].x , predicted[j].y , observations[i].x, observations[i].y);

            if(dist_to_obs < min_dist)
            {
                min_dist = dist_to_obs;

                //set the observed measurement id to the closest predictd landmark id
                observations[i].id = predicted[j].id;
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


    for (int i =0; i < num_particles; i++)
    {

        double particle_x = particles[i].x;
        double particle_y = particles[i].y;
        double particle_theta = particles[i].theta;

        //This is a vector of the observed measurements transformed to the map's coordinate system
        std::vector<LandmarkObs> trans_observations;

        //This is a vector of the landmarks that are within the sensor range only
        std::vector<LandmarkObs> predicted_landmarks;



        //Apply the transformation from vehicle to map coordinate system to the observed measurements
        for(int j=0; j<observations.size(); j++)
        {
            LandmarkObs trans_observation;
            trans_observation.id = observations[j].id;
            trans_observation.x = observations[j].x * cos(particle_theta) - observations[j].y * sin(particle_theta) + particle_x;
            trans_observation.y = observations[j].x * sin(particle_theta) + observations[j].y * cos(particle_theta) + particle_y;

            trans_observations.push_back(trans_observation);
        }

        //keep only the landmarks that are within the sensor range
        for (int j=0; j<map_landmarks.landmark_list.size(); j++)
        {
            LandmarkObs predicted_landmark;

            predicted_landmark.x = map_landmarks.landmark_list[j].x_f;
            predicted_landmark.y = map_landmarks.landmark_list[j].y_f;
            predicted_landmark.id = map_landmarks.landmark_list[j].id_i;

            if(dist(predicted_landmark.x, predicted_landmark.y, particle_x, particle_y) <= sensor_range)
            {
                predicted_landmarks.push_back(predicted_landmark);
            }

        }

        //Perform data association between the transformed observations and the map landmarks
        dataAssociation(predicted_landmarks, trans_observations);

        //use the multivariate gaussian distribution to update the particle weight
        particles[i].weight = 1.0;
        for(int j =0; j<trans_observations.size(); j++)
        {
            double predicted_x, predicted_y, observed_x, observed_y, std_x, std_y;

            observed_x = trans_observations[j].x;
            observed_y = trans_observations[j].y;
            std_x = std_landmark[0];
            std_y = std_landmark[1];


            //get the predicted landmark that is associated with the observed landmark
            for(int k =0; k<predicted_landmarks.size(); k++)
            {
                if (predicted_landmarks[k].id == trans_observations[j].id)
                {
                    predicted_x = predicted_landmarks[k].x;
                    predicted_y = predicted_landmarks[k].y;
                }
            }

            //Update the particle weights by multiplying the gaussian multivariate distributions
            particles[i].weight *= 1.0/(2* M_PI * std_x * std_y) *
                    exp(-((pow(observed_x - predicted_x, 2) / (2 * pow(std_x,2))) +
                        (pow(observed_y - predicted_y, 2) / (2 * pow(std_y,2)))));

        }

    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    default_random_engine gen;

    //This vector contains the rsampled particles
    std::vector<Particle> resampled_particles;

    //choose random index
    uniform_int_distribution<int> int_dist(0,num_particles-1);

    int index = int_dist(gen);

    double beta = 0.0;

    std::vector<double> weights;

    //collect all particles weights in the weights vector
    for (int i=0; i<num_particles; i++)
    {
        weights.push_back(particles[i].weight);
    }

    // extract the maximum weight
    double max_weight = *max_element(weights.begin(), weights.end());

    uniform_real_distribution<double> double_dist(0,max_weight);


    //run the resampling wheel
    for (int i = 0; i< num_particles; i++)
    {

        beta = beta + 2.0 * double_dist(gen);

        while (weights[index] < beta)
        {
            beta = beta - weights[index];
            index = (index+1) % num_particles;
        }
        resampled_particles.push_back(particles[index]);
    }

    particles = resampled_particles;

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
