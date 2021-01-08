#pragma once
#include <iostream>
#include <random>
#include <math.h>
#include "Naive.h" // include the Data, Params, and Distributions classes here

template<int n>
void MCMC(Prior<n> prior, Likelihood<n> lik, Proposal<n> prop, Data data, int iter, std::default_random_engine generator) {
	// Typename Data and Params are has to be included in the ".h" file that also contain
	// the Prior, Likelihood, Proposal definitions

	// Initialise params 
	Params<n> param = prior.Sample(generator);
	Params<n> param_; // empty initialisation
	double logprior = prior.LogPdf(param);
	double loglik = lik.LogPdf(data, param);
	double loglik_, logprior_, logprop, logprop_, alpha, tobeat;
	std::uniform_real_distribution<double> dis(0., 1.);

	// Iterate
	std::cout << "iter,logprior,loglik,logprop,alpha,tobeat,accept_idx,";
	std::cout << "param,param_\n"; 
	for (int i = 0; i < iter; i++) {
		param_ = prop.Sample(param, generator);

		logprior_ = prior.LogPdf(param_);
		loglik_ = lik.LogPdf(data, param_);
		logprop_ = prop.LogPdf(param_, param); // param_ -> param
		logprop = prop.LogPdf(param, param_); // param -> param_

		alpha = std::min((logprior_ - logprior) + (loglik_ - loglik) + (logprop_ - logprop), 0.);
		tobeat = MyLog(dis(generator)); // can remove tobeat and param_ in non-debug mode

		std::cout << i << "," << logprior_ << "," << loglik_ << "," << logprop << "," << alpha << "," << tobeat << "," << (alpha > tobeat);
		std::cout << "," << param << "," << param_ << "\n"; 

		if (alpha > tobeat) {
			param = param_;
			logprior = logprior_;
			loglik = loglik_;
		}
	}
}