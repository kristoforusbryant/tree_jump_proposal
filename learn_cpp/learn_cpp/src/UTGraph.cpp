#include <iostream>
#include <cassert>

template <int n>
class UTGraph {
private:
	bool adjM[n * (n - 1) / 2];
public:
	int GetSize() const { return n; }
	int GetEdgesNum() {
		int m = 0;
		for (int i = 0; i < (n * (n - 1) / 2); i++)
			m += adjM[i];
		return m;
	}
	int GetArrayIndex(int i, int j) {
			assert(i >= 0);
			assert(i <= n);
			assert(i < j); // on the right of diagonal 

			return (i * n - i * (i + 1) / 2) + j - (i + 1);
		}
	void PrintGraph() {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (j <= i) {
						std::cout << ".   ";
					}
					else {
						std::cout << adjM[GetArrayIndex(i, j)] << "   ";
					}
				}
				std::cout << std::endl;
			}
		}
	UTGraph() {
			for (int i = 0; i < n * (n - 1) / 2; i++)
				adjM[i] = false;
		}// empty constructor

	void add_edge(int i, int j) {
			adjM[GetArrayIndex(i, j)] = true;
		}
	void remove_edge(int i, int j) {
			adjM[GetArrayIndex(i, j)] = false;
		}

	UTGraph subgraph(bool* edges) {
		UTGraph<n> G;
		for (int i = 0; i < n; i++) {
			for (int j = 0; i < j; j++)
				if (edges[GetArrayIndex(i, j)]) {
					G.add_edge(i, j);
				}
		}
		return G;
	}

	void operator+(UTGraph& other) {
		for (int i = 0; i < n * (n - 1) / 2; i++)
			adjM[i] += other[i];
	}


	// other functions:  
	// from matrix to edge array, from edge list to edge array 
	// functions with UTGraph as input 
	// complete graph, complement, connected components

};