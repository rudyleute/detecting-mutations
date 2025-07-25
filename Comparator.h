#ifndef COMPARATOR_H
#define COMPARATOR_H

#include <string>

#include "Structures.h"

using namespace std;

class Comparator {
public:
	static CompRes compareMaps(
		MutationsVCF& map1,
		Mutations& map2,
		std::map<size_t, NucleoCounter> nonErrors,
		const size_t &from,
		const size_t &to
	);
};


#endif
