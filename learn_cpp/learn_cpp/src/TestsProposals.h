#pragma once

// NAIVE PROPOSAL TEST
/*
	std::default_random_engine generator;
	Prior<5> prior;
	Proposal<5> prop;
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(5,5);
	Eigen::MatrixXd U = Eigen::MatrixXd::Identity(5, 5);
	Likelihood<5> lik(3, &D); // input is address to matrix D
	Data data(6, &U);
	for (int i = 0; i < 100; ++i) {
		Params<5> p = prior.Sample(generator);
		Params<5> p_ = prop.Sample(p, generator);
		std::cout << "p: " << std::endl;
		std::cout << p << std::endl;
		std::cout << "p_: " << std::endl;
		std::cout << p_ << std::endl;
		std::cout << "PRIOR(p): " << std::endl;
		std::cout << prior.LogPdf(p) << std::endl;
		std::cout << "PRIOR(p_): " << std::endl;
		std::cout << prior.LogPdf(p_) << std::endl;
		std::cout << "LIK(p): " << std::endl;
		std::cout << lik.LogPdf(data, p) << std::endl;
		std::cout << "LIK(p_): " << std::endl;
		std::cout << lik.LogPdf(data, p_) << std::endl;
		std::cout << "PROP(p->p_): " << std::endl;
		std::cout <<  prop.LogPdf(p,p_) << std::endl;
		std::cout << "PROP(p_->p): " << std::endl;
		std::cout << prop.LogPdf(p_,p) << std::endl;

		std::cin.get();
	}

	Prior<5> prior;
	Proposal<5> prop;
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(5,5);
	Eigen::MatrixXd U = Eigen::MatrixXd::Identity(5, 5) * 5;
	Likelihood<5> lik(3, &D); // input is address to matrix D
	Data data(20, &U);

	std::default_random_engine generator;
	MCMC(prior, lik, prop, data, 7000, generator, "res_naive.csv");

*/

// CYCLE BASIS 
/* Test Basis 
	Basis<5> basis;
	for (int i = 0; i < (5 - 1) * (5 - 2) / 2; ++i) {
		basis.glist[i].PrintGraph();
		basis.glist[i].PrintEdges();
		std::cout << std::endl;

*/

/*
std::default_random_engine generator;
	Basis<5> basis;
	for (int i = 0; i < (5 - 1) * (5 - 2) / 2; ++i) {
		std::cout << "i: " << i << std::endl;
		basis.glist[i].PrintGraph();
		basis.glist[i].PrintEdges();
		std::cout << basis.glist[i] << std::endl;
	}
	Prior<5> prior(basis);
	Proposal<5> prop(basis);
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(5, 5);
	Eigen::MatrixXd U = Eigen::MatrixXd::Identity(5, 5);
	Likelihood<5> lik(3, &D); // input is address to matrix D
	Data data(6, &U);
	for (int i = 0; i < 100; ++i) {
		Params<5> p = prior.Sample(generator);
		Params<5> p_ = prop.Sample(p, generator);
		std::cout << "p: " << std::endl;
		std::cout << p << std::endl;
		std::cout << "p.Basis_idx : " << std::endl;
		p.PrintBasis();
		std::cout << "p_: " << std::endl;
		std::cout << p_ << std::endl;
		std::cout << "p_.Basis_idx : " << std::endl;
		p_.PrintBasis();
		std::cout << "PRIOR(p): " << std::endl;
		std::cout << prior.LogPdf(p) << std::endl;
		std::cout << "PRIOR(p_): " << std::endl;
		std::cout << prior.LogPdf(p_) << std::endl;
		std::cout << "LIK(p): " << std::endl;
		std::cout << lik.LogPdf(data, p) << std::endl;
		std::cout << "LIK(p_): " << std::endl;
		std::cout << lik.LogPdf(data, p_) << std::endl;
		std::cout << "PROP(p->p_): " << std::endl;
		std::cout << prop.LogPdf(p, p_) << std::endl;
		std::cout << "PROP(p_->p): " << std::endl;
		std::cout << prop.LogPdf(p_, p) << std::endl;

		std::cin.get();
	}
*/