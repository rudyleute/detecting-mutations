#include "Comparator.h"

CompRes Comparator::compareMaps(
	MutationsVCF& map1,
	Mutations& map2,
	std::map<size_t, NucleoCounter> nonErrors,
	const size_t& from,
	const size_t& to
) {
	MutationErrors diffInVCF;
	MutationErrors diffInCust;
	InBothEr errors;

	for (auto iter = map1.begin(); iter != map1.end();) {
		if (iter->first < from || iter->first >= to) {
			++iter;
			continue;
		}

		for (auto vecIter = iter->second.begin(); vecIter != iter->second.end();) {
			auto it = map2.find(iter->first);

			// If there is no index in the current implementation, report it as an error
			if (it == map2.end()) {
				diffInVCF.emplace_back(iter->first, std::get<0>(*vecIter), std::get<1>(*vecIter),
				                       nonErrors[iter->first]);
				++vecIter;
				continue;
			}

			bool found = false;

			// Iterate over the custom implementation mutations for the position iter->first
			for (auto custVecIter = it->second.begin(); custVecIter != it->second.end();) {
				// If we found the mutation in both vcf and custom implementation, remove it from both
				if (std::get<0>(*custVecIter) == std::get<0>(*vecIter) &&
					std::get<1>(*custVecIter) == std::get<1>(*vecIter)) {
					custVecIter = it->second.erase(custVecIter); // Fix: update iterator
					vecIter = iter->second.erase(vecIter);
					found = true;
					break;
				} else if (std::get<1>(*custVecIter) == std::get<1>(*vecIter)) {
					// Match only for the action - it's an error
					errors.emplace_back(iter->first, std::get<0>(*vecIter), std::get<1>(*vecIter),
					                    std::get<0>(*custVecIter), std::get<1>(*custVecIter),
					                    std::get<2>(*custVecIter));

					custVecIter = it->second.erase(custVecIter); // Fix: update iterator
					vecIter = iter->second.erase(vecIter);
					found = true;
					break;
				} else {
					++custVecIter;
				}
			}

			// If we haven't found any matching values
			if (!found) {
				if (!it->second.empty()) {
					// Safety check
					errors.emplace_back(iter->first, std::get<0>(*vecIter), std::get<1>(*vecIter),
					                    std::get<0>(it->second[0]), std::get<1>(it->second[0]),
					                    std::get<2>(it->second[0]));
				}
				vecIter = iter->second.erase(vecIter);
			}

			// Clean up empty entries in map2
			if (it->second.empty()) {
				map2.erase(it);
			}
		}

		// Clean up empty entries in map1
		if (iter->second.empty()) {
			iter = map1.erase(iter);
		} else {
			++iter;
		}
	}

	for (const auto& [key, val] : map2) {
		if (map1.find(key) == map1.end()) {
			for (const auto& val2 : val) {
				diffInCust.emplace_back(key, std::get<0>(val2), std::get<1>(val2), std::get<2>(val2));
			}
		}
	}

	return CompRes(diffInVCF, diffInCust, errors);
}
