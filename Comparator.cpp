#include "Comparator.h"

CompRes Comparator::compareMaps(
	const MutationsVCF& map1,
	const Mutations& map2,
	std::map<size_t, NucleoCounter> nonErrors,
	const size_t& from,
	const size_t& to
) {
	MutationErrors diffInVCF;
	MutationErrors diffInCust;
	InBothEr errors;

	for (const auto& [key, val1] : map1) {
		if (key < from || key >= to) continue;

		auto it = map2.find(key);
		if (it == map2.end()) diffInVCF.emplace_back(key, std::get<0>(val1), std::get<1>(val1), nonErrors[key]);
		else if (std::get<0>(it->second) != std::get<0>(val1) || std::get<1>(it->second) != std::get<1>(val1))
			errors.emplace_back(key, std::get<0>(val1), std::get<1>(val1), std::get<0>(it->second),
			                    std::get<1>(it->second), std::get<2>(it->second));
	}

	for (const auto& [key, val2] : map2) {
		if (map1.find(key) == map1.end()) {
			diffInCust.emplace_back(key, std::get<0>(val2), std::get<1>(val2), std::get<2>(val2));
		}
	}

	return CompRes(diffInVCF, diffInCust, errors);
}
