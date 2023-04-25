#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <mpi.h>

using namespace std;

// Constants
const int NUM_CITIES = 10;
const int MAX_ITERATIONS = 10000;
const int NUM_TEMPS = 100;
const double INIT_TEMP = 1000.0;
const double COOLING_RATE = 0.5;

// Define a city
struct City {
    double x, y;
};

// Generate a random TSP instance
vector<City> generate_instance(int num_cities) {
    vector<City> cities(num_cities);
    default_random_engine rng(chrono::system_clock::now().time_since_epoch().count());
    uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < num_cities; i++) {
        cities[i].x = dist(rng);
        cities[i].y = dist(rng);
    }
    return cities;
}

// Compute the Euclidean distance between two cities
double distance(const City& a, const City& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx*dx + dy*dy);
}

// Compute the length of a TSP tour
double tour_length(const vector<int>& tour, const vector<City>& cities) {
    double length = 0.0;
    for (int i = 0; i < tour.size()-1; i++) {
        length += distance(cities[tour[i]], cities[tour[i+1]]);
    }
    length += distance(cities[tour.back()], cities[tour.front()]);
    return length;
}

// Generate a random initial tour
vector<int> generate_initial_tour(int num_cities) {
    vector<int> tour(num_cities);
    for (int i = 0; i < num_cities; i++) {
        tour[i] = i;
    }
    random_shuffle(tour.begin(), tour.end());
    return tour;
}

// Simulated annealing algorithm for local search
void simulated_annealing(const vector<City>& cities, vector<int>& tour, double& temp) {
    default_random_engine rng(chrono::system_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> dist(0, NUM_CITIES-1);
    while (temp > 0.0) {
        for (int i = 0; i < MAX_ITERATIONS; i++) {
            // Generate two random indices
            int j = dist(rng);
            int k = dist(rng);
            // Swap the cities at those indices
            double old_length = tour_length(tour, cities);
            swap(tour[j], tour[k]);
            // Compute the new tour length
            // std::cout << "Iteration " << i << std::endl;
            double new_length = tour_length(tour, cities);
            // Compute the change in length
            double delta = new_length - old_length;
            // Accept the new tour with probability e^(-delta/temp)
            if (delta < 0.0 || exp(-delta/temp) > dist(rng)/(double)RAND_MAX) {
                // Accept the new tour
                break;
            } else {
                // Revert to the previous tour
                swap(tour[j], tour[k]);
            }
        }
        // Reduce the temperature
        temp *= COOLING_RATE;
    }
}

void print_cities(const vector<City>& cities) {
    for (const auto& city : cities) {
        cout << " (" << city.x << "," << city.y << ")" << endl;
    }
}

int main(int argc, char** argv) {
    // Initialize MPI
    int rank, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    // Generate the TSP instance
    vector<City> cities = generate_instance(NUM_CITIES);

    // Generate an initial tour for each process
    vector<int> tour = generate_initial_tour(NUM_CITIES);
    
    // Compute the length of the initial tour
    double length = tour_length(tour, cities);
    double best_length = length;
    
    // Initialize the temperature
    double temp = INIT_TEMP;
    
    // Perform simulated annealing with parallelization
    for (int i = 0; i < NUM_TEMPS; i++) {
        // std::cout << "Iteration " << i << " on process " << rank << std::endl;
        // Perform local search on this process
        simulated_annealing(cities, tour, temp);
        // Compute the length of the new tour
        length = tour_length(tour, cities);
        // Find the best tour among all processes
        MPI_Allreduce(&length, &best_length, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        // If this process has the best tour, broadcast it to all processes
        if (length == best_length) {
            MPI_Bcast(&tour[0], NUM_CITIES, MPI_INT, rank, MPI_COMM_WORLD);
            // std::cout << "Best tour length: " << best_length << std::endl;
        } else {
            MPI_Bcast(&tour[0], NUM_CITIES, MPI_INT, 0, MPI_COMM_WORLD);
            // std::cout << "Best tour length: " << best_length << std::endl;
        }

        // Update the temperature
        temp *= COOLING_RATE;
    }
    
    // Print the best tour length on process 0
    if (rank == 0) {
        cout << "Best tour length: " << best_length << endl;
        cout << "Cities:" << endl;
        print_cities(cities);
    }


    
    // Finalize MPI
    MPI_Finalize();
    return 0;
}

