#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "SamReader.h"

#include <iostream>
#include <string>

using namespace BamTools;

int SamReader::getAlignments(const std::string& fileName) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		return 1;
	}

	static BamAlignment alignment;
	while (reader.GetNextAlignment(alignment)) {
		std::cout << alignment.Name << std::endl;
	}
}