#include <iostream>
#include <chrono>
#include <vector>
#include <tuple>
#include <random>


constexpr int N=1000000;

std::tuple<double,double> sum1(int n, std::vector<double>& xy)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i) {
        x += xy[2*i+0];
        y += xy[2*i+1];
    }
    return {x/n, y/n};
}

std::tuple<double,double> sum2(int n, std::vector<double>& xy)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i)
        x += xy[2*i+0];
    for (int i = 0; i < n; ++i)
        y += xy[2*i+1];
    return {x/n, y/n};
}

std::tuple<double,double> sum3(int n, std::vector<double>& xy)
{
    double x = 0;
    double y = 0;
    for (int i = 0; i < n; ++i)
        x += xy[i];
    for (int i = 0; i < n; ++i)
        y += xy[n+i];
    return {x/n, y/n};
}

void init_rand(std::vector<double>& x)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(0.0, 1.0);
    for (auto& xi : x)
        xi = dist(gen);
}

template<typename F>
void timeit(const std::string& s, int ntrials, F&& lambda)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (int trial = 0; trial < ntrials; ++trial)
        lambda();
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << s << ": " << (stop-start)/ntrials << std::endl;
}

int main()
{
    std::vector<double> xy(2*N);
    init_rand(xy);
    timeit("V1", 100, [&]() { return sum1(N, xy); });
    timeit("V2", 100, [&]() { return sum2(N, xy); });
    timeit("V3", 100, [&]() { return sum3(N, xy); });
}
