#pragma once
/*

char: 1 byte  # -127 - 127
short: signed (2 bytes)
int: signed (4 bytes) # around -2b to 2b
long: signed (4 bytes) 
long: signed (8 bytes)
float: (4 bytes) 
double: (8 bytes) 
long double: (16 bytes)  approx 10e-300 to 10e300; 

*/

/*
Scratches 

	// UTGraph Scratch
	UTGraph<8> G;
	std::cout << G.GetSize() << std::endl; 
	std::cout << G.GetArrayIndex(1, 2) << std::endl;
	G.PrintGraph();

	// Array Scratch
	const int one_i = 5;
	const int two_i = 10;

	int* one = new int[one_i]; // one and two are pointers into an int array of size 5
	int* two = new int[two_i];

	for (int i = 0; i < 5; i++) {
		one[i] = i + 1;
	}

	for (int i = 0; i < 10; i++) {
		two[i] = (i + 1)*(i + 1);
	}
	PrintArray(one, one_i);
	PrintArray(two, two_i);

	//Eigen Scratch
	Eigen::Matrix4d M;
	std::cout << M << std::endl;
	std::cout << M.size() << std::endl;

	// RNG Scratch
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	MyBetaDist beta_dist(1.0, 1.0);
	double foo;

	beta_dist.PrintSamplesNum<20>(generator, 10000);

	std::cout << "Beta Distribution" << std::endl;
	for (int i = 0; i < 100; i++) {
		foo = beta_dist.Sample(generator);

		std::cout << foo << std::endl;
		std::cin.get();
	}

	std::cout << "Uniform Distribution" << std::endl;
	for (int i = 0; i < 100; i++) {
		foo = distribution(generator);
		std::cout << foo << std::endl;
		std::cout << log(foo) << std::endl;
		std::cin.get();
	}

*/