#include <iostream>
#include <random>

static class MyBetaDist {
	double a, b;
public:
	MyBetaDist(double alpha, double beta) {
		a = alpha; 
		b = beta;
	}

	double Sample(std::default_random_engine& generator) {
		std::gamma_distribution<double> X(a, 1), Y(b,1); 

		double x = X(generator);
		double y = Y(generator);
		return  (x / (x + y));
	}

	template <int breaks>
	void PrintSamplesHist(std::default_random_engine generator, int n_samples) {
		std::vector<int> v(breaks);
		
		for (int i = 0; i < n_samples; i++) {
			int rounded = (Sample(generator) * 20);
			v[rounded]++;
		}

		for (int i = 0; i < breaks; i++) 
			std::cout << (double(i) / breaks) << ":" << std::string(v[i], '*') << std::endl; 
	}

	template <int breaks>
	void PrintSamplesNum(std::default_random_engine generator, int n_samples) {
		std::vector<int> v(breaks);

		for (int i = 0; i < n_samples; i++) {
			int rounded = (Sample(generator) * 20);
			v[rounded]++;
		}

		for (int i = 0; i < breaks; i++)
			std::cout << (double(i) / breaks) << ":" << v[i] << std::endl;
	}
};