#include <iostream>
#include "dnest3/Start.h"
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest3;

int main(int argc, char** argv)
{
	Data::get_instance().load("Data/dummy_data.dat");
	MTSampler<MyModel> sampler = setup_mt<MyModel>(argc, argv);
	sampler.run();
	return 0;
}

