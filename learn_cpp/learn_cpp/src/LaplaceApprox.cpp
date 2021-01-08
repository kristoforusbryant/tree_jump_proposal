#pragma once

#include <iostream>
#include <array>
#include <cassert>
#include <vector>
#include "UTGraph.h"
#include <Eigen/Dense>


Eigen::MatrixXd Hessian(const Eigen::MatrixXd& K, const std::vector<int>& V, int delta) {
	// K is mode, V is edge list 

	int n = int(V.size() / 2);
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
	Eigen::MatrixXd K_inv = K.inverse();

	for (int e0 = 0; e0 < n; ++e0) { // e0 = (i,j) 
		for (int e1 = e0; e1 < n; ++e1) { // e1 = (l,m) 
			int i = V[2 * e0];
			int j = V[2 * e0 + 1];
			int l = V[2 * e1];
			int m = V[2 * e1 + 1];

			if (i == j && l == m) {
				H(e0, e1) = K_inv(i, l) * K_inv(i, l); 
				H(e1, e0) = H(e0, e1);
			}
			else if ((i != j && l == m) || (i == j && l != m)) {
				H(e0, e1) = K_inv(i, m) * K_inv(j, l) + K_inv(i, l) * K_inv(j, m);
				H(e1, e0) = H(e0, e1);
			}
			else {
				H(e0, e1) = 2 * (K_inv(i, m) * K_inv(j, l) + K_inv(i, l) * K_inv(j, m));
				H(e1, e0) = H(e0, e1);
			}
		}
	}

	return H * (.5) * (delta - 2);
}

double maxmin(double x) { return std::max(std::min(x, 700.), -700.); }

double LogGammaP(double x, int p) {
	double pi = 3.14159265358979323846;
	double acc = (p * (p - 1) / 4 * std::log(pi)); 
	for (int i = 0; i < p; ++i) {
		acc += std::log(tgamma(x - i/2));
	}
	return acc; 
}

double Wishart_PDF(const Eigen::MatrixXd& K, int delta, const Eigen::MatrixXd& D) {
	int p = D.rows();
	double normalising = (delta + p - 1) * p / 2 * std::log(2) + LogGammaP((delta + p - 1) / 2, p) - (delta + p - 1) / 2 * std::log(D.determinant()); 
	return (double(delta) - 2) * std::log(K.determinant()) - .5 * (K.transpose() * D).trace() - normalising;
}