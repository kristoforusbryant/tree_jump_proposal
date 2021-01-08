#pragma once

/*
	// Test choose
	std::cout << "Test Choose" << std::endl;
	std::cout << LogChoose(1, 1) << std::endl;
	std::cout << LogChoose(2, 1) << std::endl;
	std::cout << LogChoose(5, 2) << std::endl;
	std::cout << LogChoose(40, 3) << std::endl;
	std::cout << LogChoose(22, 10) << std::endl;
*/

/* UTGraph TEST 
	
	// Adding and Removing Edges 
	UTGraph<4> G0;
	UTGraph<4> G1; 
	std::cout << "G0" << std::endl; 
	G0.PrintGraph();
	std::cout << "G1" << std::endl;
	G1.PrintGraph();

	std::cout << "Adding edges to G0" << std::endl;
	G0.AddEdge(0, 1);
	G0.AddEdge(0, 2);
	G0.AddEdge(0, 3);
	G0.AddEdge(1, 2);
	G0.AddEdge(1, 3);
	G0.AddEdge(2, 3);
	G0.PrintGraph(); 
	std::cout << "Removing edges to G0" << std::endl;
	G0.RemoveEdge(0, 2);
	G0.RemoveEdge(2, 3);
	G0.RemoveEdge(1, 2);
	G0.RemoveEdge(0, 3);
	G0.PrintGraph();

	std::cout << "Adding edges to G1" << std::endl;
	G1.AddEdge(0, 2);
	G1.AddEdge(0, 3);
	G1.PrintGraph();

	// '+=' operator 
	std::cout << "Adding: G0 + G1" << std::endl;
	G0 += G1;  
	std::cout << "G0" << std::endl;
	G0.PrintGraph();
	std::cout << "G1" << std::endl;
	G1.PrintGraph();

	// Testing operator[] overload 
	std::cout << "G[0,1]: " << G0.GetArrayIndex(0,1) << std::endl;
	std::cout << "G[0,1]: " << G0(0, 1) << std::endl;
	std::cout << "G[1,2]: " << G0.GetArrayIndex(1, 2) << std::endl;
	std::cout << "G[1,2]: " << G0(1, 2) << std::endl;
	
	// Testing subgraph 
	std::array<int, 3> nodes0 = { 1,2,3 }; 
	UTGraph<3> S0; 
	G0.Subgraph<nodes0.size()>(nodes0, S0); 
	std::cout << "G0: " << std::endl; 
	G0.PrintGraph();
	std::cout << "S0: " << std::endl; 
	S0.PrintGraph(); 

	std::array<int, 2> nodes2 = { 0, 3};
	UTGraph<2> S2;
	G0.Subgraph<nodes2.size()>(nodes2, S2);
	std::cout << "G0: " << std::endl;
	G0.PrintGraph();
	std::cout << "S2: " << std::endl;
	S2.PrintGraph();

	// Neighbours
	std::vector<int> nbh;
	nbh = G0.Neighbours(3);
	std::cout << "Neighbours of 3: " << std::endl;
	for (int i = 0; i < int(nbh.size()); i++)
		std::cout << nbh[i] << std::endl;
	std::cout << std::endl;

	nbh = G0.Neighbours(0);
	std::cout << "Neighbours of 0: " << std::endl;
	for (int i = 0; i < int(nbh.size()); i++)
		std::cout << nbh[i] << std::endl;
	std::cout << std::endl;

	// Complete Graph
	UTGraph<4> G1 = CompleteGraph<4>();
	UTGraph<8> G2 = CompleteGraph<8>();
	std::cout << "G1: " << std::endl;
	G1.PrintGraph();
	std::cout << "G2: " << std::endl;
	G2.PrintGraph();

	// Complement Graph
	G1 = ComplementGraph<4>(G0);
	std::cout << "G0: " << std::endl;
	G0.PrintGraph();
	std::cout << "G1: " << std::endl;
	G1.PrintGraph();

	UTGraph<8> RG0;
	UTGraph<15> RG1;

	std::default_random_engine generator;

	for (int i = 0; i < 100; i++) {
		// Random
		RG0 = RandomGraph<8>(generator);
		//RG1 = RandomGraph<15>(generator);
		std::cout << "RG0" << std::endl;
		RG0.PrintGraph();
		RG0.PrintEdges();
		//std::cout << "RG1" << std::endl;
		//RG1.PrintGraph();

		// DFUtils
		std::array<bool, 8> visited = { false };
		DFUtils<8>(RG0, 0, visited);
		PrintArray<std::array<bool, 8>>(visited, 8);

		// ConnectedComps
		std::array<int, 8> concomp = { 0 };
		ConnectedComponents<8>(RG0, concomp);
		PrintArray<std::array<int, 8>>(concomp, 8);

		// PrimeComps 
		MySet<8> P(true);
		MySet<8> R; 
		MySet<8> X;
		std::vector<MySet<8>> primecomp; 
		PrimeComponents<8>(RG0, P, R, X, primecomp);
		for (int i = 0; i < int(primecomp.size()); i++)
			primecomp[i].Print();

		std::cin.get();
	}

*/

/* MySet
	std::cout << std::endl << "TESTING MYSET..." << std::endl;
	MySet<8> S;
	S.Add(1);
	S.Add(6);
	std::cout << "S = " <<  S << std::endl;
	MySet<8> S_ = S.AsArray();
	std::cout << "S_ =" << S_ << std::endl;
	S_.Complement();
	std::cout << "S_ =" << S_ << std::endl;

	S_ = S.GetComplement();
	std::cout << "S_ =" << S_ << std::endl;

	std::vector<int> v = S.AsIntVector();
	std::cout << "{ ";
	for (size_t i = 0; i < v.size(); i++) {
		std::cout << v[i] << ", ";
	}
	std::cout << "}" << std::endl;

*/


/*
	// Test Prior
	Prior beta_prior_ev(5, 5);
	std::cout << "Beta(5,5)" << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(.5)) << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(0)) << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(2)) << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(-1)) << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(.7)) << std::endl;
	std::cout << beta_prior_ev.LogPdf(Params(.3)) << std::endl;
	std::cout << beta_prior_ev.LogRatio(Params(.3), Params(.7)) << std::endl;

	// Test Likelihood
	Likelihood likelihood_ev(20);
	std::cout << "Binom(20,p)" << std::endl;
	std::cout << likelihood_ev.LogPdf(Data(10), Params(.5)) << std::endl;
	std::cout << likelihood_ev.LogPdf(Data(10), Params(.7)) << std::endl;
	std::cout << likelihood_ev.LogPdf(Data(10), Params(.3)) << std::endl;
	std::cout << likelihood_ev.LogRatio(Data(10), Params(.3), Params(.7)) << std::endl;

	std::cout << likelihood_ev.LogPdf(Data(5), Params(.5)) << std::endl;
	std::cout << likelihood_ev.LogPdf(Data(5), Params(.7)) << std::endl;
	std::cout << likelihood_ev.LogPdf(Data(5), Params(.3)) << std::endl;
	std::cout << likelihood_ev.LogRatio(Data(5), Params(.3), Params(.7)) << std::endl;

	// Test Prior, Likelihood, Proposal
	Prior prior(1,2);
	Likelihood lik(5);
	Proposal prop(.5);

	std::default_random_engine generator;
	// Prior
	for (int i = 0; i < 100; i++) {
		Params p = prior.Sample(generator);
		std::cout << "(" << p.GetTheta() << ", " << exp(prior.LogPdf(p)) <<"); ";
	}
	std::cout << std::endl;

	// Proposal
	Params p(.4);
	Data d(3);
	for (int i = 0; i < 100; i++) {
		Params p_ = prop.Sample(p, generator);
		std::cout << "(" << p_.GetTheta() << ", Prior: " << exp(prior.LogPdf(p_)) << ", Lik: " << exp(lik.LogPdf(d, p_)) << ", Prop: " << exp(prop.LogPdf(p, p_)) << "); " << std::endl;
	}
*/

/* MCMC TEST
	Data data(2);
	Prior prior(1, 1);
	Likelihood lik(5);
	Proposal prop(.7);

	std::default_random_engine generator;
	MCMC(prior, lik, prop, data, 5000, generator);
*/

/* PLAY AROUND WITH EIGEN 
 	Eigen::MatrixXd m0 = Eigen::MatrixXd::Random(3, 3);
	Eigen::MatrixXd m1 = Eigen::MatrixXd::Random(3, 3);

	std::cout << "m0 =" << std::endl << m0 << std::endl;
	std::cout << "m1 =" << std::endl << m1 << std::endl;

	std::cout << std::endl  << "Updating Matrices ..." << std::endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			m0(i, j) = static_cast<int>(10.0 * abs(m0(i, j)));
			m1(i, j) = static_cast<int>(10.0 * abs(m1(i, j)));
		}
	}

	std::cout << std::endl << "MATRIX MULT ..." << std::endl;
	std::cout << "m0 =" << std::endl << m0  << std::endl;
	std::cout << "m1 =" << std::endl << m1 << std::endl;
	std::cout << "m0 * m1 =" << std::endl << m0 * m1 << std::endl;

	std::cout << std::endl << "CHECK FOR BROADCASTING ..." << std::endl;
	std::cout << "m0 =" << std::endl << m0 / 2 << std::endl;

	std::cout << std::endl  << "INVERSES ..." << std::endl;
	m0 = Eigen::MatrixXd::Identity(3, 3) * 5;
	m1 += Eigen::MatrixXd::Identity(3, 3) * 20; 
	Eigen::MatrixXd m2 = Eigen::MatrixXd::Identity(3, 3); 
	m2(0, 0) = 0.0; 

	std::cout << "m0 =" << std::endl << m0 << std::endl;
	std::cout << "m0_inv =" << std::endl << m0.inverse() << std::endl;
	std::cout << "m1 =" << std::endl << m1 << std::endl;
	std::cout << "m1_inv =" << std::endl << m1.inverse() << std::endl;
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << "m1(" << i << "," << j << ") = " << m1(i, j) << std::endl;
		}
	}

	std::cout << "m2 =" << std::endl << m2 << std::endl;
	std::cout << "THIS SUPPOSED TO FAIL" << std::endl;
	std::cout << "m2_inv =" << std::endl << m2.inverse() << std::endl;

	std::cout << std::endl << "CHECK THAT MATMUL WITH THE WRONG DIMENSIONS FAIL..." << std::endl;
	Eigen::MatrixXi r0 = Eigen::MatrixXi::Random(2, 4);
	Eigen::MatrixXi r1 = Eigen::MatrixXi::Random(4, 4);

	std::cout << "r0 =" << std::endl << r0 << std::endl;
	std::cout << "r1 =" << std::endl << r1 << std::endl;
	std::cout << "r0 * r1 =" << std::endl << r0 * r1 << std::endl;
	// std::cout << "r1 * r0 =" << std::endl << r1 * r0 << std::endl; // this gives assertion runtime error

	std::cout << std::endl << "CHECK INDEXING ..." << std::endl;
	std::vector<int> indices0 = { 0,2 };
	std::vector<int> indices1 = { 1,3 };

	std::cout << "r1({0,2}, {1,3}) = " << std::endl << r1(indices0, indices1) << std::endl;
	r1(indices0, Eigen::all) = r0;
	std::cout << "r1= " << std::endl << r1 << std::endl;

*/

/* Test Multivariate Normal Distribution Generator
	
	/*int dim = 4;
	int N = 10;
	Eigen::MatrixXd covar(dim, dim);
	covar << 1, .0, .0, .0,
		.0, 1, .0, .0,
		.0, .0, 1, .0,
		.0, .0, .0, 1;

	normal_random_variable sample{ covar };
	Eigen::VectorXd s;
	Eigen::MatrixXd dat(N, dim);

	for (int i = 0; i < N; i++) {
		s = sample();
		dat({ i }, Eigen::all) = s;
		for (int j = 0; j < dim; j++) {
			std::cout << s[j] << ", ";
		}
		std::cout << std::endl;
	}

	std::cout << "### Data Matrix D ###" << std::endl; 
	std::cout << dat << std::endl;

*/

// LAPLACE APPROX

/*  Mode
	UTGraph<4> G1; // Empty graph
	int delta = 5;
	Eigen::MatrixXd D(4, 4);
	D << 1, .0, .0, .0,
		.0, 1, .0, .0,
		.0, .0, 1, .0,
		.0, .0, .0, 1;

	std::cout << "Mode: " << std::endl;
	Eigen::MatrixXd K = Mode(G1, delta, D);
	std::cout << K << std::endl;

	double mode_pdf = Wishart_PDF(K, delta, D); 
	
	Eigen::MatrixXd RM0;
	Eigen::MatrixXd R;
	std::default_random_engine generator;

	for (int i = 0; i < 1000; ++i) {

		// Random Diagonal Matrix
		RM0 = Eigen::MatrixXd::Random(4, 4).cwiseAbs();
		std::cout << "RM0: " << std::endl;
		std::cout << RM0 << std::endl;

		R = RM0.diagonal().asDiagonal(); 
		std::cout << "R: " << std::endl;
		std::cout << R << std::endl;
		
		
		std::cout << "Wishart Sample PDF: " << Wishart_PDF(R, delta, D) << std::endl;
		std::cout << "Wishart Mode PDF: " << mode_pdf << std::endl;
	
		std::cin.get(); 
	}
*/


/* ConstrainedConc

	UTGraph<8> RG0;
	Eigen::MatrixXd RM0;
	Eigen::MatrixXd L;
	Eigen::MatrixXd R;
	std::default_random_engine generator;

	for (int i = 0; i < 100; i++) {
		// Random Graph
		RG0 = RandomGraph<8>(generator);
		std::cout << "RG0: " << std::endl;
		RG0.PrintGraph();
		RG0.PrintEdges();

		// Random Symmetric Positive Definite Matrix
		RM0 = Eigen::MatrixXd::Random(8, 8).cwiseAbs();
		L = RM0 * RM0.transpose();
		std::cout << "L: " << std::endl;
		std::cout << L << std::endl;


		// Constrained Conc
		R = ConstrainedConc(RG0, L);
		std::cout << "R: " << std::endl;
		std::cout << R << std::endl;
		std::cout << "R.inverse(): " << std::endl;
		std::cout << R.inverse() << std::endl << std::endl;

		std::cin.get();
	}

*/

/* Hessian 
	UTGraph<5> RG0;
	Eigen::MatrixXd RM0;
	Eigen::MatrixXd K;
	Eigen::MatrixXd D;
	Eigen::MatrixXd H;
	


	std::default_random_engine generator;

	for (int i = 0; i < 100; i++) {
		// reset edges 
		std::vector<int> edges;
		for (int i = 0; i < 5; i++) { // fill the diagonals 
			edges.push_back(i);
			edges.push_back(i);
		}

		// Random Graph
		RG0 = RandomGraph<5>(generator);
		std::cout << "RG0: " << std::endl;
		RG0.PrintGraph();
		RG0.PrintEdges();
		RG0.EdgeList(edges); 
		std::cout << "edges size: " << int(edges.size() / 2) << std::endl; 
		std::cout << "edges: ";
		for (int i = 0; i < int(edges.size() / 2); ++i) {
			std::cout << "(" << edges[2 * i] << "," << edges[2 * i + 1] << "), ";
		}
		std::cout << std::endl;

		
		// Random Symmetric Positive Definite Matrices
		RM0 = Eigen::MatrixXd::Random(5, 5).cwiseAbs();
		K = RM0 * RM0.transpose();
		std::cout << "K_inv: " << std::endl;
		std::cout << K.inverse() << std::endl;

		// Constrained Conc 
		H = Hessian(K, edges, 10);
		std::cout << "H: " << std::endl;
		std::cout << H << std::endl;

		std::cin.get();
	}

*/

/* Laplace Approx Altogether 
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
			std::cout << prop.LogPdf(p, p_) << std::endl;
			std::cout << "PROP(p_->p): " << std::endl;
			std::cout << prop.LogPdf(p_, p) << std::endl;

			std::cin.get();
		}


	Params<5> p;
	Eigen::MatrixXd D = Eigen::MatrixXd::Identity(5, 5);
	Eigen::MatrixXd U = Eigen::MatrixXd::Identity(5, 5) * 5;
	Likelihood<5> lik(3, &D); // input is address to matrix D
	Data data(20, &U);
	Eigen::MatrixXd D_star;

	D_star = ConstrainedConc<5>(p.GetG(), D + *(data.GetU()));
	std::cout << "p: " << std::endl;
	std::cout << p << std::endl;
	std::cout << "LIK(p): " << std::endl;
	std::cout << lik.LogPdf(data, p) << std::endl;

	std::cout << "D: " << std::endl;
	std::cout << D << std::endl;
	std::cout << LaplaceApprox(p.GetG(), 3, D) << std::endl;
	std::cout << "D_star: " << std::endl;
	std::cout << D_star << std::endl;
	std::cout << LaplaceApprox(p.GetG(), 100, D_star) << std::endl;

	for (int i = 0; i < 10; ++i) {
		p.FlipEdge(i);
		D_star = ConstrainedConc<5>(p.GetG(), D + *(data.GetU()));
		std::cout << "p: " << std::endl;
		std::cout << p << std::endl;
		std::cout << "LIK(p): " << std::endl;
		std::cout << lik.LogPdf(data, p) << std::endl;
	}

*/
