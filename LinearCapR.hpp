#pragma once

#include <vector>
#include <string>

#include "types.hpp"

using namespace std;

class LinearCapR{
public:
	LinearCapR(int beam_size) : beam_size(beam_size){}


private:
	const int beam_size;
};