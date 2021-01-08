#include <iostream>
#include <cassert>
#include <random>
#include <math.h>
#include <Eigen/Dense>
#include "UTGraph.h"
#include "LaplaceApprox.h"
#include "MultiNorm.h"
#include "MCMC.h"


int main() {
	Prior<11> prior;
	Proposal<11> prop;
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(11,11);
	Eigen::MatrixXd U = Eigen::MatrixXd::Identity(11,11) * 5;
	Likelihood<11> lik(3, &D); // input is address to matrix D
	Data data(20, &U);
	
	std::default_random_engine generator;
	MCMC(prior, lik, prop, data, 100, generator);
	
	

	/*UTGraph<4> RG0 = RandomGraph<4>(generator);
	std::cout << "RG0: " << std::endl;
	RG0.PrintGraph();
	RG0.PrintEdges();
	bool* foo  = RG0.GetArray();*/
	

	//std::cin.get();
}