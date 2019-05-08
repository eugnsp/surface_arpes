#include "simulator.hpp"

#include <mkl_service.h>

#include <cstdlib>
#include <exception>
#include <iostream>

int main(int argc, const char** argv)
{
	try
	{
		std::cout << "1D Poisson solver\n"
				  << "=================\n"
				  << std::endl;

		Simulator sim;
		sim.run(argc, argv);
		::mkl_free_buffers();
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception!\n" << e.what() << '\n';
		return EXIT_FAILURE;
	}
	catch (...)
	{
		std::cout << "Exception!\n";
		return EXIT_FAILURE;
	}

	std::cout << "Simulation completed.\n";
	return EXIT_SUCCESS;
}
