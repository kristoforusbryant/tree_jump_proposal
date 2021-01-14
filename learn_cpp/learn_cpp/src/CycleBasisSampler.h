#pragma once
#include <iostream>
#include <cmath> // same as math.h 
#include <vector>
#include <random>
#include "UTGraph.h"
double MyLog(double n);
double LogChoose(int n, int k);

/*
SUPPORT TO SAMPLED GRAPHS FROM CYCLE BASIS 
*/

// Define Basis 
template<int n>
struct CycleBasis {
	UTGraph<n> T;
	UTGraph<n> glist[(n - 1) * (n - 2) / 2];

	CycleBasis() {
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
class CycleParams {
	UTGraph<n> G;
	std::array<bool, (n - 1)* (n - 2) / 2> Basis_idx{ false };
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

	friend std::ostream& operator<<(std::ostream& out, const CycleParams& param) {
		out << param.G;
		return out;
	}

	template<int n>
	friend class CyclePrior;
};

// Uniform Prior on the Space of Graphs 
template<int n>
class CyclePrior {
	CycleBasis<n> basis;
public:
	CyclePrior(CycleBasis<n>& basis_) {
		basis = basis_;
	}

	// Sample 
	CycleParams<n> Sample(std::default_random_engine& generator) {
		CycleParams<n> p;
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
	double LogPdf(const CycleParams<n>& param) {
		return -((n - 1) * (n - 2) / 2) * MyLog(2);
	}
};
