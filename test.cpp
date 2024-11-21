#include <iostream>
#include <iterator>
#include <random>

int main() {
    // Example data
    std::vector<double> data = {1., 2., 3., 4., 5., 6.};

    // Define random generator with Gaussian distribution
    const double mean = 0.0;//均值
    const double stddev = 0.1;//标准差
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (auto& x : data) {
        x = x + dist(generator);
    }

    // Output the result, for demonstration purposes
    std::copy(begin(data), end(data), std::ostream_iterator<double>(std::cout, " "));
    std::cout << "\n";

    return 0;
}