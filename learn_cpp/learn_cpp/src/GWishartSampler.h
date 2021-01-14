#pragma once
#include <iostream>
#include <array>
#include <cassert>
#include <vector>
#include <Eigen/Dense>
#include "UTGraph.h"
#include "Wishart.h"
#include "MultiNorm.h"

template<int n>
class GWishartLikelihoodSampler {
private: 
	Eigen::MatrixXd K; 
	UTGraph<n> G; 
	int delta; 
	Eigen::MatrixXd D_inv; 
	int N; 
	std::vector<MySet<n>> primecomps;
	double* w;
	
public: 
	GWishartLikelihoodSampler(const UTGraph<n>& G_, int delta_, const Eigen::MatrixXd& D_, int N_ = 50) {
		K = Eigen::MatrixXd::Identity(n, n);
		G = G_; 
		delta = delta_; 
		D_inv = D_.inverse(); 
		N = N_; 

		// Generating PrimeComps from G
		MySet<n> P(true);
		MySet<n> R;
		MySet<n> X;
		{
			UTGraph<n> Gtemp = G; // Copy contents of G, destruct after using PrimsComponents
			PrimeComponents<n>(Gtemp, P, R, X, primecomps); // returns a vector of MySet (array of bools)
		}
		
		/*for (size_t i = 0; i < primecomps.size(); i++) {
			std::cout << primecomps[i] << ", ";
		}
		std::cout << std::endl;*/

		for (int i = 0; i < N; i++) StepOnce();

	}

	void StepOnce() {
		// temporary variables
		Eigen::MatrixXd B_;
		Eigen::MatrixXd C_;
		Eigen::MatrixXd D_;
		for (size_t j = 0; j < primecomps.size(); j++) {
			std::vector<int> temp = primecomps[j].AsIntVector();
			std::vector<int> temp_ = MySet<n>(primecomps[j].GetComplement()).AsIntVector();
			int s = temp.size();

			Eigen::VectorXi c = Eigen::Map<Eigen::VectorXi>(temp.data(), s);
			Eigen::VectorXi c_ = Eigen::Map<Eigen::VectorXi>(temp_.data(), s);

			B_ = K(c, c_); // rows in the clique, columns not  // GetComplement already returns std::array<bool,n>
			C_ = K(c_, c); // rows not in the clique, columns is
			D_ = K(c_, c_); // rows and columns not in the clique

			double* sigma = new double[s*s];
			for (int i = 0; i < s; ++i) {
				for (int j = 0; j < s; ++j) {
					sigma[i + j * s] = D_inv(temp[i], temp[j]);
				}
			}

			w = wishart_sample(s, delta, sigma);
			Eigen::MatrixXd W = Eigen::MatrixXd::Zero(s, s);

			for (int i = 0; i < s; ++i) {
				for (int j = i; j < s; ++j) {
					W(i, j) = w[i + s * j];
					W(j, i) = w[i + s * j];
				}
			}

			K(c, c) = W.inverse() + B_ * D_.inverse() * C_;

			free(sigma);
		}
	}

	Eigen::VectorXd operator()(std::default_random_engine& generator) {
		StepOnce();
		normal_random_variable sample{ K.inverse() };
		return sample(generator);
	}
};
