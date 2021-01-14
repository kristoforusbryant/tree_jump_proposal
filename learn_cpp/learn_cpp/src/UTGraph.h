#pragma once

#include <iostream>
#include <cassert>
#include <array>
#include <random>

struct Edges {
	std::array<int, 2> arr;
	Edges(int i, int j) {
		arr[0] = i;
		arr[1] = j;
	}
	void PrintEdges() {
		std::cout << "(" << arr[0] << "," << arr[1] << ")" << std::endl;
	}
};

template <int n>
class UTGraph {
private:
	static_assert(n > 1, "Number of vertices < 2");
	bool adjM[n * (n - 1) / 2];
public:
	int GetSize() const { return n; }
	bool GetArrayValue(int i) const { return adjM[i]; }
	void SetArrayValue(int i, bool b) { adjM[i] = b;  }

	int GetEdgesNum() const {
		int m = 0;
		for (int i = 0; i < (n * (n - 1) / 2); i++)
			m += adjM[i];
		return m;
	}
	int GetArrayIndex(int i, int j) const {
		assert(i >= 0);
		assert(i <= n);
		assert(i < j); // on the right of diagonal 

		return (i * n - i * (i + 1) / 2) + j - (i + 1);
	}
	bool GetElement(int i, int j) const {
		assert(i != j);
		int min = std::min(i, j); 
		int max = std::max(i, j);
		return adjM[GetArrayIndex(min, max)];
	}
	void PrintGraph() const {
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
	void PrintEdges() const {
		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				if (GetElement(i,j)) std::cout << "(" << i << "," << j << ") "; 
			}
		}
		std::cout << std::endl;
	}
	UTGraph() {
		for (int i = 0; i < n * (n - 1) / 2; i++)
			adjM[i] = false;
	} // empty constructor
	

	UTGraph(const UTGraph<n>& other) {
		for (int i = 0; i < n * (n - 1) / 2; i++)
			adjM[i] = other.GetArrayValue(i);
	} // copy constructor

	void AddEdge(int i, int j) {
		adjM[GetArrayIndex(std::min(i,j), std::max(i,j))] = true;
	}
	void RemoveEdge(int i, int j) {
		adjM[GetArrayIndex(std::min(i, j), std::max(i, j))] = false;
	}

	void BinaryAddToArray(int i) { 
		adjM[i] = adjM[i] ^ 1;
	}

	void AddCircle() {
		for (int i = 0; i < n * (n - 1) / 2; i += (n - 1 - i)) { adjM[i] = true; }
	}

	template <int n>
	void Neighbours(int idx, std::array<bool, n> &nb) {
		assert(idx < n);
		for (int i = 0;  i < n; i++) {
			if (idx == i) continue;
			if (GetElement(idx, i)) nb[i] = true;
		}
	}

	template<int n>
	bool FindCycle(int now, int prev, std::vector<int>& visited) { // find the first cycle via DFS starting from the node "now"
		visited.push_back(now);
		std::array<bool, n> nb{};
		Neighbours(now, nb); 
		/*std::cout << "now: " << now << std::endl; 
		std::cout << "neighbours: [";
		for (int i = 0; i < n; ++i) {
			std::cout << nb[i] << ", " ; 
		}
		std::cout << "]" << std::endl;;*/
		
		for (int i = 0; i < n; ++i) { // a neighbour, but not previous
			if (nb[i] && (i != prev)) {
				if (std::find(visited.begin(), visited.end(), i) != visited.end()) {
					visited.push_back(i);
					std::reverse(visited.begin(), visited.end());
					while ( visited.back() != i) {
						visited.pop_back();
					};
					return true;
				}
				else { 
					if (FindCycle<n>(i, now, visited)) { return true; };
				}
			}
		}
		visited.pop_back(); 
		return false;
	}


	template <int n_nodes>
	void Subgraph(std::array<int, n_nodes> arr, UTGraph<n_nodes> & other) { // not tested
		for (int i = 0; i < n_nodes; i++) {
			for (int j = i+1; j < n_nodes; j++) {
				std::cout << "(i,j): (" <<  i << "," << j << ")" << std::endl;
				std::cout << "arr(i,j): (" << arr[i] << "," << arr[j] << ")" << std::endl;
				if (GetElement(arr[i], arr[j])) other.AddEdge(i,j); 
			}
		}
	}

	void EdgeList(std::vector<int>& edges) const { // make sure this is how we should pass vector by reference
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (GetElement(i, j)) {
					edges.push_back(i);
					edges.push_back(j);
				}
			}
		}
	}
	
	int operator()(int i, int j) const {
		return GetArrayIndex(i, j);
	}
	void operator+=(const UTGraph& other) { 
		for (int i = 0; i < n * (n - 1) / 2; i++)
			adjM[i] = adjM[i] | other.adjM[i];
	} 

	void BinaryAddition(const UTGraph& other) {
		for (int i = 0; i < n * (n - 1) / 2; i++)
			adjM[i] = adjM[i] ^ other.adjM[i];
	}


	friend std::ostream& operator<<(std::ostream& out, const UTGraph& G) {
		out << "[ ";
		for (int i = 0; i < int(n * (n - 1) / 2); i++)
			out << G.adjM[i] << " ";
		out << "]";
		return out;
	}

};



template <int n>
UTGraph<n> CompleteGraph() {
	UTGraph<n> G;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			G.AddEdge(i, j); 
		}
	}
	return G; 
}

template <int n>
UTGraph<n> ComplementGraph(const UTGraph<n>& other) {
	UTGraph<n> G = CompleteGraph<n>();
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (other.GetElement(i, j)) G.RemoveEdge(i, j); 
		}
	}
	return G; 
}

template <int n>
void DFUtils(const UTGraph<n>& G, int x, std::array<bool, n>& visited) { // replace visited in place
	visited[x] = true;
	std::array<bool, n> nb = { false };
	G.Neighbours<n>(x, nb);
	for (int i = 0; i < n; i++) {
		if (nb[i] & !visited[i]) {
			DFUtils(G, i, visited);
		}
	}
}



template <int n>
void ConnectedComponents(const UTGraph<n>& G, const std::array<int, n>& comps) {
	int comps_i = 0;
	std::array<bool, n> visited = { false }; 
	std::array<bool, n> visited_temp = { false };
	for (int i = 0; i < n; i++) {
		if (visited[i]) { continue; }
		visited_temp.fill(false);
		DFUtils<n>(G, i, visited_temp);

		for (int j = 0; j < n; j++) {
			if (visited_temp[j]) { comps[j] = comps_i; }
			visited[j] = visited[j] | visited_temp[j];
		}
		comps_i++;
	}
}

// Random Graph  
template <int n> 
UTGraph<n> RandomGraph(std::default_random_engine &generator) {
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	UTGraph<n> G; 
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (distribution(generator) < .5) G.AddEdge(i, j);
		}
	}
	std::cout << std::endl;
	return G; 
}


// other functions:  
	// from matrix to edge array, from edge list to edge array 