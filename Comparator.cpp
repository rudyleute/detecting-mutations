#include "Comparator.h"

compRes Comparator::compareMaps(const Mutations& map1,
					const Mutations& map2) {
	std::vector<std::tuple<size_t, char, char>> difference1;
	std::vector<std::tuple<size_t, char, char>> difference2;
	std::vector<std::tuple<size_t, char, char, char, char>> errors;

	for (const auto& [key, val1] : map1) {
		auto it = map2.find(key);
		if (it == map2.end()) difference1.emplace_back(key, val1.first, val1.second);
		else if (it->second != val1)
			errors.emplace_back(key, val1.first, val1.second, it->second.first,
								it->second.second);
	}

	for (const auto& [key, val2] : map2) {
		if (map1.find(key) == map1.end()) difference2.emplace_back(key, val2.first, val2.second);
	}

	return compRes(difference1, difference2, errors);
}
