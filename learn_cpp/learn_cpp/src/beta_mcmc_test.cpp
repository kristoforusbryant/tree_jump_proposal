#include <iostream>
#include <cmath> // same as math.h 
#include <vector>
#include <random>
double MyLog(double n);
double LogChoose(int n, int k);

// Define Parameter Type 
class Params {
	double theta; 
	
public: 
	Params(double t) {
		theta = t; 
	}
	void PrintParams() {
		std::cout << theta << std::endl; 
	}
	double GetTheta() {
		return theta;
	}

	friend std::ostream& operator<<(std::ostream& out, const Params & param) {
		out << param.theta;
		return out;
	}
};

// Define Data Type 
class Data {
	int n_success; 
public: 
	Data(int s){
		n_success = s;
	}
	void PrintData() {
		std::cout << n_success << std::endl;
	}
	int GetS() {
		return n_success;
	}
};


// Beta Prior Pdf 
// Is it better to represent distributions as functions or objects??? 
class Prior{
	double a, b, log_beta_f;
public: 
	Prior(double alpha, double beta) {
		a = alpha;
		b = beta;
		log_beta_f =(MyLog(std::tgamma(a)) + MyLog(std::tgamma(b)) - MyLog(std::tgamma(a + b)));
	}

	// Sample 
	Params Sample(std::default_random_engine& generator) {
		std::gamma_distribution<double> X(a, 1), Y(b, 1);

		double x = X(generator);
		double y = Y(generator);
		Params p = Params(x / (x + y));
		return  p;
	}

	// Evaluate
	double LogPdf(Params p) {
		if ((p.GetTheta() > 1) || (p.GetTheta() < 0)) { return -700; }
		else { return (a - 1.) * MyLog(p.GetTheta()) + (b - 1.) * MyLog(1. - p.GetTheta()) - log_beta_f; }
	} 
};

// Binomial Likelihood for 1 Observation 
class Likelihood {
	int n; 
public: 
	Likelihood (int num) {
		n = num; 
	}

	double LogPdf(Data data, Params p){
		if ((p.GetTheta() > 1) || (p.GetTheta() < 0)) { return -700; }
		else {
			return LogChoose(n, data.GetS()) + double(data.GetS()) * MyLog(p.GetTheta()) +
				(double(n) - double(data.GetS())) * MyLog(double(1) - double(p.GetTheta()));
		}
	}
};


// Uniform Proposal of length "a" around the previous parameter
class Proposal {
	double a; 
public: 
	Proposal(double a_) {
		a = a_;
	}

	// Sample 
	Params Sample(Params p, std::default_random_engine& generator) {
		std::uniform_real_distribution<double> dis(-a, +a);

		Params p_ = Params(dis(generator) + p.GetTheta());
		return  p_;
	}

	// Evaluate
	double LogPdf(Params p0, Params p1) { // i.e. p0 -> p1 
		if (p1.GetTheta() > p0.GetTheta() + a || p1.GetTheta() < p0.GetTheta() - a) { return -700; }
		else { return 0 - MyLog(2 * a); }
	}

};