#include <cmath>
#include <cassert>

double MyLog(double n) {
	if (log(n) > 700) { return 700; }
	if (log(n) < -700) { return -700; }
	return log(n);
}

double LogChoose(int n, int k) {
	assert(k >= 0); 
	if (k < 1) { return 0;}
	if (n < k) { return -700;}
	
	double prod = 0;
	for (int i = 0; i < k; i++)
		prod += MyLog(double(n) - double(i)) - MyLog(double(k) - double(i));

	return prod;
}

// Log factorial