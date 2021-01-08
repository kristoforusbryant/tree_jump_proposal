#pragma once
#include <iostream>
#include <array>
#include <cassert>
#include <vector>
#include "UTGraph.h"
#include <Eigen/Dense>


template <int n>
class MySet {
private:
	std::array<bool, n> set;
public:
	MySet() {
		set = { false };
	}
	MySet(bool b) {
		for (int i = 0; i < n; i++)
			set[i] = b;
	}
	MySet(std::array<bool, n> a) { // copy constructor from array 
		for (int i = 0; i < n; i++)
			set[i] = a[i];
	}
	MySet(const MySet& other) { // copy constructor from another MySet
		for (int i = 0; i < n; i++)
			set[i] = other.set[i];
	}

	void Print() const {
		std::cout << "{ "; 
		for (int i = 0; i < n; i++)
			std::cout << set[i] << ", ";
		std::cout << "}"<< std::endl;
	}

	int GetSize() const { return n; }
	int GetCount() const {
		int count = 0;
		for (int i = 0; i < n; i++) {
			if (set[i]) count++; 
		}
		return count;
	}
	bool IsEmpty() const {
		for (int i = 0; i < n; i++) {
			if (set[i]) { return false; }
		}
		return true;
	}
	int Pop() {
		assert(!IsEmpty());
		for (int i = 0; i < n; i++) {
			if (set[i]) {
				set[i] = false;
				return i;
			}
		}
	}
	void Add(int i) { 
		assert(0 <= i && i < n);
		set[i] = true; 
	}
	void operator+=(MySet& other) {
		for (int i = 0; i < n; i++)
			set[i] = set[i] | other.set[i];
	}// union 
	void operator-=(MySet& other) {
		for (int i = 0; i < n; i++)
			set[i] = set[i] & !(other.set[i]);
	}// setminus
	void Intersect(MySet& other) {
		for (int i = 0; i < n; i++)
			set[i] = set[i] & other.set[i];
	}
	void Intersect(std::array<bool, n>& other) {
		for (int i = 0; i < n; i++)
			set[i] = set[i] & other[i];
	}
	
	// setcomplement 
	void Complement() {
		for (int i = 0; i < n; i++) {
			set[i] = (set[i] ^ 1);
		}
	}

	std::array<bool, n> GetComplement() const {
		std::array<bool, n> new_set;
		for (int i = 0; i < n; i++) {
			new_set[i] = set[i] ^ true;
		}
		return new_set; 
	}

	std::array<bool, n> AsArray() const {
		return set; 
	}

	std::vector<int> AsIntVector() const {
		std::vector<int> v;
		for (int i = 0; i < n; i++) {
			if (set[i]) v.push_back(i);
		}
		return v;
	}

	friend std::ostream& operator<<(std::ostream& out, const MySet& myset) {
		out << "{ ";
		for (int i = 0; i < n; i++) {
			out << myset.set[i] << ", ";
		}
		out << "}";
		return out;
	}
};

template<int n>
void PrimeComponents(UTGraph<n> &G, MySet<n>& P, MySet<n>& R, MySet<n>& X, std::vector<MySet<n>> &primes) { // Bron Kerbosch  
	if (P.IsEmpty() & X.IsEmpty()) {
		primes.push_back(R);
	}

	while (!P.IsEmpty()) {
		int v = P.Pop();
		// Is fixed array implementation the most space efficient?
		std::array<bool, n> nb = { false }; 
		G.Neighbours<n>(v, nb);

		MySet<n> P_(P); // using copy constructor
		P_.Intersect(nb); // in-place operation
		MySet<n> R_(R);
		R_.Add(v); 
		MySet<n> X_(X);
		X_.Intersect(nb);

		PrimeComponents(G, P_, R_, X_, primes); 
		X.Add(v);
	}
}

template<int n>
Eigen::MatrixXd Mode(const UTGraph<n>& G, int delta, const Eigen::MatrixXd& D, int N = 50) {
	Eigen::MatrixXd L = D / (delta - 2);  
	Eigen::MatrixXd K = Eigen::MatrixXd::Identity(n, n); 
	Eigen::MatrixXd K_; // a copy to allow computation of differences between iterations 

	// Generating PrimeComps from G
	MySet<n> P(true);
	MySet<n> R;
	MySet<n> X;
	std::vector<MySet<n>> primecomps;
	{
		UTGraph<n> G_ = G; // Copy contents of G, destruct after using PrimsComponents
		PrimeComponents<n>(G_, P, R, X, primecomps); // returns a vector of MySet (array of bools)
	}
	
	/*for (size_t i = 0; i < primecomps.size(); i++) {
		std::cout << primecomps[i] << ", ";
	}
	std::cout << std::endl;*/

	// temporary variables
	Eigen::MatrixXd B_; 
	Eigen::MatrixXd C_;
	Eigen::MatrixXd D_;
	for (int i = 0; i < N; i++) {
		K_ = K;
		for (size_t j = 0; j < primecomps.size(); j++) {
			std::vector<int> temp = primecomps[j].AsIntVector(); 
			std::vector<int> temp_ = MySet<n>(primecomps[j].GetComplement()).AsIntVector();
			
			Eigen::VectorXi c = Eigen::Map<Eigen::VectorXi>(temp.data(), temp.size()); 
			Eigen::VectorXi c_ = Eigen::Map<Eigen::VectorXi>(temp_.data(), temp_.size());

			B_ = K(c, c_); // rows in the clique, columns not  // GetComplement already returns std::array<bool,n>
			C_ = K(c_, c); // rows not in the clique, columns is
			D_ = K(c_, c_); // rows and columns not in the clique
			
			K(c, c) = L(c, c).inverse() + B_ * D_.inverse() * C_;
		}
		if ((K_ - K).cwiseAbs().maxCoeff() < 1e-50) break; // max norm
	}
	// TODO: Give runtime warning if the value does not converge 
	return K;
}

template<int n>
Eigen::MatrixXd ConstrainedConc(const UTGraph<n>& G, const Eigen::MatrixXd& L, Eigen::MatrixXd M = Eigen::MatrixXd::Identity(n,n), int N = 50) {
	// The second cyclic algorithm in (Speed and Kiiveri, 1986) 
	// G is the Graph 
	// L is the matrix that matches existent edges of G
	// M is the matrix that matches missing edges in G
	
	Eigen::MatrixXd K = M.inverse();
	Eigen::MatrixXd K_; // a copy to allow computation of differences between iterations 

	// Generating PrimeComps from G
	MySet<n> P(true);
	MySet<n> R;
	MySet<n> X;
	std::vector<MySet<n>> primecomps;
	// TODO: consider doing this once and passing it as reference to both ConstrainedConc and Mode
	{
		UTGraph<n> G_ = G; // Copy contents of G, destruct after using PrimsComponents
		PrimeComponents<n>(G_, P, R, X, primecomps); // returns a vector of MySet (array of bools)
	}

	
	for (int i = 0; i < N; i++) {
		K_ = K;
		for (size_t j = 0; j < primecomps.size(); j++) {
			std::vector<int> temp = primecomps[j].AsIntVector();
			std::vector<int> temp_ = MySet<n>(primecomps[j].GetComplement()).AsIntVector();

			Eigen::VectorXi c = Eigen::Map<Eigen::VectorXi>(temp.data(), temp.size());
			Eigen::VectorXi c_ = Eigen::Map<Eigen::VectorXi>(temp_.data(), temp_.size());

			Eigen::MatrixXd B_c = L(c, c);

			Eigen::MatrixXd R_c_inv = K(c, c).inverse();
			Eigen::MatrixXd R_cc_ = K(c, c_); 
			Eigen::MatrixXd R_c_c = K(c_, c);
			Eigen::MatrixXd R_c_c_ = K(c_, c_);

			K(c, c) = B_c;
			K(c_, c) = R_c_c * R_c_inv * B_c;
			K(c, c_) = B_c * R_c_inv * R_cc_;
			K(c_, c_) = R_c_c_ - R_c_c * R_c_inv * (Eigen::MatrixXd::Identity(c.rows(), c.rows()) - B_c * R_c_inv) * R_cc_;
		}
		if ((K_ - K).cwiseAbs().maxCoeff() < 1e-50) break; // max norm
	}
	// TODO: Give runtime warning if the value does not converge 
	return K;
}

Eigen::MatrixXd Hessian(const Eigen::MatrixXd& K, const std::vector<int>& V, int delta);

double maxmin(double x);

template<int n>
double LaplaceApprox(const UTGraph<n>& G, int delta, const Eigen::MatrixXd& D) {
	// Laplace Approximation of the G-Wishart as outlined by (Lenkoski and Dobra, 2011) 
	Eigen::MatrixXd K = Mode(G, delta, D); 
	std::vector<int> edges;
	for (int i = 0; i < n; i++) { // fill the diagonals 
		edges.push_back(i); 
		edges.push_back(i);
	}
	G.EdgeList(edges); 
	/*std::cout << "G: " << std::endl; 
	std::cout << G << std::endl;
	std::cout << "Mode: " << std::endl;
	std::cout << K << std::endl;*/
	double h = -.5 * ((K.transpose() * D).trace() - (delta - 2) * std::log(K.determinant())); 
	/*std::cout << h << std::endl; 
	std::cout << (edges.size() / 2 / 2) * std::log(2 * 3.14159265358979323846) << std::endl;
	std::cout << -(.5) * std::log(Hessian(K, edges, delta).determinant()) << std::endl;
	std::cout << Hessian(K, edges, delta) << std::endl;*/

	return h + (edges.size() / 2 / 2) * std::log(2 * 3.14159265358979323846) - (.5) * std::log(Hessian(K, edges, delta).determinant());
}


double LogGammaP(double x, int p); 
double Wishart_PDF(Eigen::MatrixXd& K, int delta, Eigen::MatrixXd& D);

// TODO: Support for evaluating G-Wishart exactly for decomposable graphs

