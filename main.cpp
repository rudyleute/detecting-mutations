#include <filesystem>
#include <iostream>

#include "SamReader.h"

int main(int argc, char *argv[]) {
  const std::string fullPath = std::filesystem::current_path().parent_path().string() + '/' + argv[1];
  auto alignments = SamReader::getAlignments(fullPath);
}