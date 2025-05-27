#ifndef COMPARATOR_H
#define COMPARATOR_H

#include <string>

#include "Structures.h"

using namespace std;

class Comparator {
public:
	static CompRes compareMaps(
		const MutationsVCF& map1,
		const Mutations& map2,
		const size_t &from,
		const size_t &to
	);
};


#endif
