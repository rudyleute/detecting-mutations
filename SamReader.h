#ifndef SAMREADER_H
#define SAMREADER_H
#include <string>

class SamReader {
public:
	static int getAlignments(const std::string& fileName);
};

#endif //SAMREADER_H
