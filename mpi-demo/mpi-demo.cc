#include <mpi.h>

#include <print>
#include <vector>
#include <span>
#include <random>

/**
 * ## Hello world
 *
 * This one is relatively straightforward: each rank says hello by
 * printing to the standard output.  Note that the output lines can
 * appear in whatever order they want!
 *
 */
void hello(int size, int rank)
{
    std::println("Hello from {} / {}", rank, size);
}


/**
 * ## Send, recv, and sendrecv
 *
 * The `MPI_Send` operation is blocking, which means that under some
 * circumstances (like those we described in class), we can get into a
 * deadlock if there is a cycle of processors trying to send to each
 * other who mutually can't make it to their receive operations until
 * the sends complete.  *But* whether this actually results in deadlock
 * or not depends on the size of the messages and the amount of buffering
 * that the MPI implementation (and operating system) provides!  Here
 * we test this out, running with incrementally larger messages until
 * we actually see a deadlock condition.
 *
 */
void deadlock_attempt(int size, int rank)
{
    if (size != 2) {
        std::println("Should be run with two ranks");
        return;
    }

    int num_ints = 1;
    while (1) {
        MPI_Status status;
        std::vector<int> sendbuf(num_ints);
        std::vector<int> recvbuf(num_ints);
        MPI_Send(sendbuf.data(), num_ints, MPI_INT, 1-rank,
                 0, MPI_COMM_WORLD);
        MPI_Recv(recvbuf.data(), num_ints, MPI_INT, 1-rank,
                 0, MPI_COMM_WORLD, &status);
        std::println("{}: Received from {} ({})",
                     rank, status.MPI_SOURCE, num_ints);
        num_ints *= 10;
    }
}


/**
 * The pairing of an overlapping send and receive is pretty common,
 * enough to deserve special support with `MPI_Sendrecv`.  Though
 * the overall send-receive pair blocks, the sub-operations (the
 * individual sends and recieves) run concurrently, so that we won't
 * deadlock even if there is a communication cycle.  We test this
 * out here.
 *
 */
void deadlock_nonattempt(int size, int rank)
{
    if (size != 2) {
        std::println("Should be run with two ranks");
        return;
    }

    int num_ints = 1;
    for (int k = 0; k < 7; ++k) {
        MPI_Status status;
        std::vector<int> sendbuf(num_ints);
        std::vector<int> recvbuf(num_ints);
        MPI_Sendrecv(sendbuf.data(), num_ints, MPI_INT, 1-rank, 0,
                     recvbuf.data(), num_ints, MPI_INT, 1-rank, 0,
                     MPI_COMM_WORLD, &status);
        std::println("{}: Received from {} ({})", 
                     rank, status.MPI_SOURCE, num_ints);
        num_ints *= 10;
    }
}


/**
 * ## Computing all-to-all interactions
 *
 * The typical all-to-all interaction looks like
 * $$
 *   f(x_i) = \sum_{j} k(x_i, x_j) c_j
 * $$
 * where $k$ is an interaction kernel.  For some kernels, there is
 * a faster way to compute this than explicitly evaluating all the
 * sums, but that is not our interest here.  We want to show off
 * some MPI calls by doing things in a brute-force way.
 *
 * The "pass it around" protocol discussed in the slides manages
 * memory scalability by saying that each rank is assigned "ownership"
 * of disjoint pieces of the array of points $x_i$.  A rank is responsible
 * for accumulating $f(x_i)$ for the points that it owns; but to do so,
 * it needs to see all the other points.  We manage this with a
 * "pass the dish" protocol: we keep a set of particle data that we
 * compute interactions against at each phase, and rotate those interactions
 * around the ring.
 *
 */

/**
 * ### The main data structure
 *
 * Let's set up a class for all the communication buffers and
 * indexing structures that we need.
 *
 */
class A2ASim {
public:
    A2ASim(int npoints);
    void compute();

private:
    // -- Number of procs and local rank in MPI communicator
    int nproc;
    int rank;

    // -- Offsets and counts of points owned by each rank
    std::vector<int> offsets; // Start offsets for each rank in global idx
    std::vector<int> counts;  // Counts of data owned at each rank
    int max_count;            // Max counts owned by any rank

    // -- Accumulators and points (local)
    std::vector<double> results; // Accumulators for local results
    std::vector<double> points;  // Local points

    // -- Buffers for pass-it-around protocol
    std::vector<double> rpoints; // Remote points for current computation
    std::vector<double> sendbuf; // Buffer for sending data
    std::vector<double> recvbuf; // Buffer for receiving data
    MPI_Request req[2];          // Send and receive requests

    // Helper functions
    void generate_points();
    void print_results();
    void start_sendrecv(int phase);
    void end_sendrecv();
    void add_interact(std::span<double> other_pts);
};


// Initialize an a2a_sim_t data structure
A2ASim::A2ASim(int npoints)
{
    // Get size and rank info from communicator
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Set up offsets array
    offsets.resize(nproc+1);
    counts.resize(nproc);
    for (int i = 0; i <= nproc; ++i)
        offsets[i] = (i * npoints)/nproc;

    // Compute counts from offsets
    max_count = 0;
    for (int i = 0; i < nproc; ++i) {
        counts[i] = offsets[i+1]-offsets[i];
        if (counts[i] > max_count)
            max_count = counts[i];
    }

    // Set up rest of data
    results.resize(counts[rank]);
    points.resize(counts[rank]);
    rpoints.resize(max_count);
    sendbuf.resize(max_count);
    recvbuf.resize(max_count);
}


/**
 * ### Scattering and gathering
 *
 * In order to avoid the complexities of parallel random number generation,
 * we are going to create all the initial points at rank 0 and then *scatter*
 * them to the other processors.
 *
 */
void A2ASim::generate_points()
{
    if (rank == 0) {
        // Initialize the global points array
        int ntotal = offsets[nproc];
        std::vector<double> all_points(ntotal);
        std::print("Initialize array:");
        for (int i = 0; i < ntotal; ++i) {
            all_points[i] = drand48();
            std::print(" {:g}", all_points[i]);
        }
        std::println("");
        
        // Distribute points to all ranks (scatter)
        MPI_Scatterv(all_points.data(),
                     counts.data(), offsets.data(), MPI_DOUBLE,
                     points.data(), counts[rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);
    } else {
        // Distribute points to all ranks (scatter)
        MPI_Scatterv(nullptr,
                     counts.data(), offsets.data(), MPI_DOUBLE,
                     points.data(), counts[rank], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);
    }
    
    // Print what we got at each rank
    std::print("On rank {} ({}-{}): ", rank,
                 offsets[rank], offsets[rank+1]);
    for (int i = 0; i < counts[rank]; ++i)
        std::print(" {:g}", points[i]);
    std::println("");
}


/**
 * Conversely, at the end we will *gather* all the results to rank 0
 * for output.
 *
 */
void A2ASim::print_results()
{
    if (rank == 0) {
        
        // Allocate space to receive results
        int ntotal = offsets[nproc];
        std::vector<double> all_results(ntotal);

        // Gather data from all ranks to print out
        MPI_Gatherv(results.data(), counts[rank], MPI_DOUBLE,
                    all_results.data(), counts.data(),
                    offsets.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // Print results at rank 0
        std::print("Results: ");
        for (auto result : all_results)
            std::print(" {:g}", result);
        std::println("\n");
        
    } else {
        MPI_Gatherv(results.data(), counts[rank], MPI_DOUBLE,
                    nullptr, counts.data(),
                    offsets.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);
    }
}


/**
 * ### Passing data
 *
 * In the main loop, we'll send from each processor to the next processor
 * mod the total number of processors.  We do this with a nonblocking
 * send/receive operation.  Our communication thus has two phases:
 *
 * 1. We copy data into the send buffer and start the send/receive pair.
 * 2. We wait on message completions and copy data out of the recv buffer.
 *
 * In between these two phases, we can do computations that do not involve
 * either of the buffers being communicated.
 *
 */
void A2ASim::start_sendrecv(int phase)
{
    // Previous and next processors in a ring
    int rprev = (rank + nproc - 1) % nproc;
    int rnext = (rank + 1) % nproc;

    // Copy current remote points into sendbuf space
    std::copy(sendbuf.begin(), sendbuf.end(), rpoints.begin());

    // Start nonblocking send/recieve
    MPI_Isend(sendbuf.data(), max_count, MPI_DOUBLE, rnext, phase,
              MPI_COMM_WORLD, req+0);
    MPI_Irecv(recvbuf.data(), max_count, MPI_DOUBLE, rprev, phase,
              MPI_COMM_WORLD, req+1);
}


void A2ASim::end_sendrecv()
{
    // Wait for requests to complete
    MPI_Waitall(2, req, MPI_STATUSES_IGNORE);

    // Copy data from recvbuf into rpoints for next phase of computation
    std::copy(rpoints.begin(), rpoints.end(), recvbuf.begin());
}

/**
 * ### Core computation
 *
 * We start with the computational kernel, with a $1/r^2$ interaction function
 * (the same type that we would see in electrostatic or gravitational
 * interactions).  We use unit weights on all the interactions.
 *
 */
void A2ASim::add_interact(std::span<double> other_pts)
{
    int nmine = counts[rank];
    for (int i = 0; i < nmine; ++i) {
        double fi = 0.0;
        double xi = points[i];
        for (auto xj : other_pts) {
            double dij = xi-xj;
            if (dij != 0)
                fi += 1.0/(dij*dij);
        }
        results[i] += fi;
    }
}


/**
 * ### The main computation
 *
 * Having wrapped all the MPI calls in previous helper functions,
 * the main routine to generate points, do all the interactions,
 * and gather and print results can now be fairly terse.
 *
 */
void A2ASim::compute()
{
    generate_points();
    start_sendrecv(0);
    for (int phase = 0; phase < nproc; ++phase) {
        int rdata = ((rank-phase) + nproc) % nproc;
        add_interact(std::span{rpoints.begin(), rpoints.begin()+counts[rdata]});
        end_sendrecv();
        if (phase < nproc-1)
            start_sendrecv(phase+1);
    }
    print_results();
}


/**
 * ## The main routine
 */
int main(int argc, char** argv)
{
    int world_size;
    int world_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Hello world exercise
    hello(world_size, world_rank);

    // Only run this if you want to deadlock!
    // deadlock_attempt(world_size, world_rank);

    // Or not
    // deadlock_nonattempt(world_size, world_rank);

    // Compute pairwise interactions
    A2ASim sim(5);
    sim.compute();

    MPI_Finalize();
    return 0;
}
