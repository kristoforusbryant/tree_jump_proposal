#pragma once
#include <iostream>
#include <cmath> // same as math.h 
#include <vector>
#include <random>
#include <Eigen/Dense>
#include "UTGraph.h"
#include "LaplaceApprox.h"
double MyLog(double n);
double LogChoose(int n, int k);

/*
PRIOR: UNIFORM
PROPOSAL: CYCLE BASIS
LIKELIHOOD: MARGINAL P(G | DATA) APPROXIMATED BY LAPLACE APPROX
*/

// Define Basis 
template<int n>
struct Basis {
	UTGraph<n> T; 
	UTGraph<n> glist[(n-1) * (n-2) / 2];

	Basis() {
		for (int i = 0; i < n - 1; ++i) {
			T.AddEdge(i, i + 1);
		}

		int idx = 0;
		for (int i = 0; i < n - 1; ++i) {
			for (int j = i + 2; j < n; ++j) {
				UTGraph<n> T_ = T;
				T_.AddEdge(i, j);

				std::vector<int> visited;
				T_.FindCycle<n>(0, -1, visited);

				for (size_t k = 0; k < (visited.size() - 1); ++k) {
					glist[idx].AddEdge(visited[k], visited[k + 1]);
				}
				idx++;
			}
		}

	}
};


// Define Parameter Type 
template<int n>
class Params {
	UTGraph<n> G;
	std::array<bool, (n - 1)* (n - 2) / 2> Basis_idx { false };
public:
	UTGraph<n> GetG() const {
		return G;
	}

	void FlipEdge(int i) { // using the array index directly 
		G.BinaryAddToArray(i);
	}
	
	void PrintGraph() const {
		G.PrintGraph();
	}

	void PrintEdges() const {
		G.PrintEdges();
	}

	void PrintBasis() const {
		std::cout << "[ ";
		for (int i = 0; i < (n - 1) * (n - 2) / 2; ++i) {
			std::cout << Basis_idx[i] << ", ";
		}
		std::cout << "]" << std::endl;
	}

	friend std::ostream& operator<<(std::ostream& out, const Params& param) {
		out << param.G;
		return out;
	}

	template<int n>
	friend class Prior;
	template<int n>
	friend class Likelihood;
	template<int n>
	friend class Proposal;
};

// Define Data Type 
class Data {
	int N;
	Eigen::MatrixXd* U; // TODO: consider changing this to reference
public:
	Data(int N_, Eigen::MatrixXd* U_) {
		N = N_;
		U = U_;
	}
	int Getn() const { return N; }
	Eigen::MatrixXd* GetU() const { return U; }

	template<int n>
	friend class Likelihood;
};


// Define MCMC Output
template <int n>
struct DataOutput {
	int iter;
	double logprior;
	double loglik;
	double logprop;
	double alpha;
	double tobeat;
	Params<n> param;
	Params<n> param_;

	DataOutput() {
		iter = -1;
		logprior = -1.;
		loglik = -1.;
		logprop = -1.;
		alpha = -1.;
		tobeat = -1.;
	}

	DataOutput(int iter, double logprior, double loglik, double logprop, double alpha, double tobeat, Params<n> param, Params<n> param_) {
		this->iter = iter;
		this->logprior = logprior;
		this->loglik = loglik;
		this->logprop = logprop;
		this->alpha = alpha;
		this->tobeat = tobeat;
		this->param = param;
		this->param_ = param_;
	}

	friend std::ostream& operator<<(std::ostream& out, const DataOutput& data) {
		out << data.iter << "," << data.logprior << "," << data.loglik << "," << data.logprop << "," << data.alpha << "," << data.tobeat << "," << (data.alpha > data.tobeat);
		out << "," << data.param << "," << data.param_ << "\n";
		return out;
	}
};



// Uniform Prior on the Space of Graphs 
template<int n>
class Prior {
	Basis<n> basis; 
public:
	Prior(Basis<n>& basis_) {
		basis = basis_;
	}

	// Sample 
	Params<n> Sample(std::default_random_engine& generator) {
		Params<n> p;
		std::bernoulli_distribution dis(.5);
		for (int i = 0; i < int((n - 1) * (n - 2) / 2); ++i) {
			if (dis(generator)) {
				p.G.BinaryAddition(basis.glist[i]);
				p.Basis_idx[i] = true;
			}
		}
		return p;
	}

	// Evaluate
	double LogPdf(const Params<n>& param) {
		return -((n - 1) * (n - 2) / 2) * MyLog(2);
	}
};

// Marginal Likelihood estimate as ratio of G-Wishart Normalisation Constants
template<int n>
class Likelihood {
	int delta;
	Eigen::MatrixXd* D; // pointer to the matrix D 
public:
	Likelihood(int delta_, Eigen::MatrixXd* D_) {
		delta = delta_;
		D = D_;
	}

	double LogPdf(const Data& data, const Params<n>& param) {
		Eigen::MatrixXd L = *D + *(data.U);
		Eigen::MatrixXd D_star = ConstrainedConc<n>(param.G, L); // either make both see the same memory address, or copy G into the this scope 
		/*std::cout << "IG(delta, D_star): " << std::endl;
		std::cout << LaplaceApprox<n>(param.GetG(), delta + data.Getn(), D_star) << std::endl;
		std::cout << "IG(delta, D): " << std::endl;
		std::cout << LaplaceApprox<n>(param.GetG(), delta, *D) << std::endl;*/

		return LaplaceApprox<n>(param.G, delta + data.N, D_star) - LaplaceApprox<n>(param.G, delta, *D);
	}
};


// Uniform Proposal
template<int n>
class Proposal {
	Basis<n> basis;
public:
	Proposal(Basis<n>& basis_) {
		basis = basis_;
	}

	// Sample 
	Params<n> Sample(const Params<n>& p, std::default_random_engine& generator) {
		Params<n> p_ = p;
		std::uniform_int_distribution<int> dis(0, ((n - 1) * (n - 2) / 2) - 1);
		int i = dis(generator);
		p_.G.BinaryAddition(basis.glist[i]); // Test this 
		p_.Basis_idx[i] = p.Basis_idx[i] ^ true;
		return p_;
	}

	// Evaluate
	double LogPdf(const Params<n>& p0, const Params<n>& p1) { // i.e. p0 -> p1 
		return 0;  // Symmetric proposal 
	}

};
