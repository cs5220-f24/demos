#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char** argv)
{
    mpi::environment env;
    mpi::communicator world;

    std::cout << "Hello from processor " << world.rank();
    std::cout << " out of " << world.size() << " processors." << std::endl;

    return 0;
}
