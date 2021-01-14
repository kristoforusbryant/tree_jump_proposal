#include <iostream>
#include <fstream>
#include <cassert>
#include <random>
#include <math.h>
#include <Eigen/Dense>
#include "UTGraph.h"
#include "LaplaceApprox.h"
#include "GWishartSampler.h"
#include "CycleBasisSampler.h"
#include "MCMC.h" 
#include "ExperimentsNaive.h"

#include<chrono>
#include<thread>

struct Timer
{
	std::chrono::time_point<std::chrono::steady_clock> start, end;
	std::chrono::duration<float> duration;

	Timer() {
		std::cout << "Starting Timer! \n";
		start = std::chrono::high_resolution_clock::now();
		;
	}

	~Timer() {
		end = std::chrono::high_resolution_clock::now();
		duration = end - start;

		float s = duration.count();
		std::cout << "Timer took " << s << "s \n";
	}

};

int main() {
	std::ios_base::sync_with_stdio(false); 
	std::cin.tie(NULL);
	
	{
		Timer time;
		experiment0<11>(7500);

	}

	/*
	{
		Timer time;
		experiment1<11>(7500);

	}

	{
		Timer time;
		experiment2<11>(7500);

	}

	{
		Timer time;
		experiment3<11>(7500);

	}*/
	//Prior<5> prior;
	//Proposal<5> prop;
	//Eigen::MatrixXd D = Eigen::MatrixXd::Identity(5, 5);
	//Eigen::MatrixXd U = Eigen::MatrixXd::Identity(5, 5) * 5;
	//Likelihood<5> lik(3, &D); // input is address to matrix D
	//Data data(20, &U);

	//std::default_random_engine generator;
	//MCMC(prior, lik, prop, data, 100, generator, "res_naive.csv");


	//std::cin.get();
}