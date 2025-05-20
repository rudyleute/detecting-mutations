#include <iostream>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "FilesReader.h"

#include <list>
#include <boost/icl/interval_map.hpp>
#include <htslib/sam.h>
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>

using namespace BamTools;
using namespace std;

char FilesReader::complementSymbol(const char& n) {
	switch (n) {
	case 'A': return 'T';
	case 'T': return 'A';
	case 'G': return 'C';
	case 'C': return 'G';
	case '-': return '-';
	default: return 'N';
	}
}

std::map<char, size_t> FilesReader::nucleoMapping = {
	{'A', 0}, {'G', 1}, {'C', 2}, {'T', 3}, {'-', 4}
};

std::map<size_t, char> FilesReader::rNucleoMapping = [] {
	std::map<size_t, char> reverse;
	for (const auto& [key, value] : FilesReader::nucleoMapping) {
		reverse[value] = key;
	}
	return reverse;
}();

std::map<std::string, std::tuple<size_t, size_t, size_t>> FilesReader::cigarIndices;

void FilesReader::reverseComplement(std::string& read) {
	reverse(read.begin(), read.end());
	transform(read.begin(), read.end(), read.begin(), complementSymbol);
}

string FilesReader::getRefGen(const string& fileName) {
	string refGen;

	ifstream fin(fileName);
	string line;
	while (getline(fin, line)) {
		if (line.empty() || line[0] == '>') continue;
		refGen += line;
	}
	fin.close();

	return refGen;
}

CigarString FilesReader::getCigarString(bam1_t *b) {
	auto cigar = bam_get_cigar(b);
	size_t n_cigar = b->core.n_cigar;

	CigarString cigarData;
	for (int i = 0; i < n_cigar; ++i) {
		int op_len = bam_cigar_oplen(cigar[i]);
		char op_char = bam_cigar_opchr(cigar[i]);
		cigarData.emplace_back(op_char, op_len);
	}

	return cigarData;
}

string FilesReader::getRead(bam1_t* b, bool isReversed) {
	auto* seq = bam_get_seq(b);
	int32_t seqLen = b->core.l_qseq;

	string read(seqLen, 'a');
	// if (!isReversed) {
	// 	for (int i = 0; i < seqLen; i++) {
	// 		switch (bam_seqi(seq, i)) {
	// 			case 1: read[i] = 'A'; break;
	// 			case 2: read[i] = 'C'; break;
	// 			case 4: read[i] = 'G'; break;
	// 			case 8: read[i] = 'T'; break;
	// 			case 15: read[i] = 'N'; break;
	// 			default: read[i] = '?'; break;
	// 		}
	// 	}
	// } else {
	// 	for (int i = 0; i < seqLen; i++) {
	// 		switch (bam_seqi(seq, i)) {
	// 			case 1: read[i] = 'T'; break;
	// 			case 2: read[i] = 'G'; break;
	// 			case 4: read[i] = 'C'; break;
	// 			case 8: read[i] = 'A'; break;
	// 			case 15: read[i] = 'N'; break;
	// 			default: read[i] = '?'; break;
	// 		}
	// 	}
	// }

	for (int i = 0; i < seqLen; i++) {
		switch (bam_seqi(seq, i)) {
			case 1: read[i] = 'A'; break;
			case 2: read[i] = 'C'; break;
			case 4: read[i] = 'G'; break;
			case 8: read[i] = 'T'; break;
			case 15: read[i] = 'N'; break;
			default: read[i] = '?'; break;
		}
	}

	return read;
}

string FilesReader::getExpandedRead(string read, CigarString &cigar) {
	string expandedRead;
	if (cigar.begin()->first == 'S' || cigar.begin()->first == 'H') {
		read = read.substr(cigar.begin()->second);
		cigar.erase(cigar.begin());
	}

	auto cigarLast = next(cigar.end(),  -1);
	if (cigarLast->first == 'S' || cigarLast->first == 'H') {
		read = read.substr(0, read.size() - cigarLast->second);
		cigar.erase(cigarLast);
	}

	size_t readFrom = 0;
	size_t curInd = 0;
	for (const auto &cigarElem: cigar) {
		if (cigarElem.first == 'D') {
			expandedRead += read.substr(readFrom, curInd - readFrom);
			readFrom = curInd;

			expandedRead += string(cigarElem.second, '-');
		} else curInd += cigarElem.second;
	}

	if (curInd != readFrom) expandedRead += read.substr(readFrom, curInd - readFrom);

	return expandedRead;
}

std::map<size_t, std::pair<char, char>> FilesReader::readReferenceCsv(const string& fileName) {
	std::ifstream file(fileName);
	if (!file.is_open()) {
		std::cerr << "Failed to open the file.\n";
		throw runtime_error("Error opening file " + fileName);
	}

	std::map<size_t, std::pair<char, char>> dataMap;
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string fStr, sStr, tStr;

		if (!std::getline(ss, fStr, ',')) continue;
		if (!std::getline(ss, sStr, ',')) continue;
		if (!std::getline(ss, tStr, ',')) continue;

		if (fStr.empty() || tStr.empty()) continue;

		//TODO remove this line
		// if (fStr[0] == 'I') continue;

		char fCol = fStr[0];
		size_t sCol = std::stoull(sStr);
		char tCol = tStr[0];

		dataMap[sCol] = std::make_pair(tCol, fCol);
	}

	return dataMap;
}

size_t FilesReader::getReferenceLength(const std::string& fileName) {
	BamReader reader;

	if (!reader.Open(fileName)) {
		std::cerr << "Error opening file " << fileName << std::endl;
		throw runtime_error("Error opening file " + fileName);
	}

	return reader.GetReferenceData()[0].RefLength;
}

AlignmentMaps FilesReader::getAlignments(const std::string& fileName, const size_t& from, const size_t& to,
                                         InsertionMap prevIterInsertions) {
	samFile *in = sam_open(fileName.c_str(), "r");
	bam_hdr_t *header = sam_hdr_read(in);

	hts_idx_t *idx = sam_index_load(in, fileName.c_str());

	hts_itr_t *iter = sam_itr_querys(idx, header, ("NC_001416:" + std::to_string(from) + "-" + std::to_string(to)).c_str());
	bam1_t *b = bam_init1();

	std::map<size_t, set<pair<string, string>>> startingPos;
	std::map<size_t, pair<NucleoCounter, set<string>>> insertions = std::move(prevIterInsertions);
	std::map<size_t, pair<NucleoCounter, set<string>>> insertionsOOB;

	while (sam_itr_next(in, iter, b) >= 0) {
		//Check whether the sequence has been aligned to the reference genome
		if (b->core.flag & BAM_FUNMAP) continue;

		size_t pos = b->core.pos;
		string name = bam_get_qname(b);
		CigarString cigarExpanded = getCigarString(b);

		cout << name << " " << (b->core.flag & BAM_FREVERSE) << endl;
		const string read = getRead(b, b->core.flag & BAM_FREVERSE);
		string expandedRead = getExpandedRead(read, cigarExpanded);

		if (name == "NC-001416-mutated_12202_aligned_166_R_24_12985_15") {
			cout << 1;
		}

		//Aligned position for the string with cut out insertions (for the direct substitution and deletion analysis)
		size_t readSDStart = pos;
		size_t startPos = pos;

		// If we analyze parts of the string that are not within the window, we can find ourselves in the situation when indices are analyzed more than once with different content
		// check whether the starting index is within the window
		// if so, iterate over the cigar string clipping it in a way that would cover only the part of the read within the window
		// store the index of the last analyzed element in the Cigar string
		// if the starting index is not within the window check the array with stored indices in order to figure out what the position that should be started with and remove the index

		size_t readOffset = 0;
		size_t cigarOffset = 0;
		if (cigarIndices.find(name) != cigarIndices.end()) {
			auto aux = cigarIndices[name];
			const auto end = next(cigarExpanded.begin(), get<0>(cigarIndices[name]));
			cigarExpanded.erase(cigarExpanded.begin(), end);
			cigarExpanded.begin()->second -= get<1>(cigarIndices[name]);
			expandedRead = expandedRead.substr(get<2>(cigarIndices[name]));
			cigarOffset = get<0>(cigarIndices[name]);
			readOffset = get<2>(cigarIndices[name]);
			readSDStart = from;
			startPos = from;

			cigarIndices.erase(name);
		}
		else if (startPos < from) continue;
		//Since the insertions are cut out, if the first non-clipped part of the string requires an insertion
		//the starting index should take that into consideration
		if (cigarExpanded.begin()->first == 'I') readSDStart += cigarExpanded.begin()->second;
		//Base index for the insertion - should take care of the potential clipping

		size_t insertionsFound = 0;
		string noInsertionsRead;
		size_t readFromIndex = 0;
		size_t curReadIndex = 0;

		for (auto iter = cigarExpanded.begin(); iter != cigarExpanded.end(); iter++) {
			const size_t refGenIndex = startPos - insertionsFound + curReadIndex;
			const bool isInsertion = iter->first == 'I';

			//if we get out of the boundaries of the defined window
			if (refGenIndex + iter->second >= to) {
				//We can cover the insertions even if they get out of the window as it does not increase the indices for S and D
				const size_t left = to - refGenIndex;
				const size_t end = isInsertion ? iter->first : left;
				auto* curMap = &insertions;

				if (isInsertion) {
					noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
					readFromIndex = curReadIndex + iter->second;
				}

				for (size_t i = 0; i != end; i++) {
					if (i == left) {
						curMap = &insertionsOOB;
					}

					if (isInsertion) {
						(*curMap)[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
						(*curMap)[refGenIndex + i].second.insert(name);
					}
					else if ((*curMap)[refGenIndex + i].second.find(name) == (*curMap)[refGenIndex + i].
					                                                                   second.
					                                                                   end()) {
						(*curMap)[refGenIndex + i].first.increase('-');
					}
				}

				if (!isInsertion) {
					//We want to store the position in the iteration string, length covered and position in the expanded string
					//in order to start the analysis of the string within the next window quicker
					cigarIndices[name] = make_tuple(distance(cigarExpanded.begin(), iter) + cigarOffset, left,
					                                          curReadIndex + left + readOffset);
					curReadIndex += left;
					break;
				}

				insertionsFound += iter->second;
				curReadIndex += iter->second;
				continue;
			}


			//TODO ensure that only the parts of the string within the window are analyzed in order to avoid the same indices to be analyzed several times and other errors
			if (isInsertion) {
				noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
				readFromIndex = curReadIndex + iter->second;
			}

			for (size_t i = 0; i != iter->second; i++) {
				if (refGenIndex + i == 1056) {
					cout << 2;
				}
				if (isInsertion) {
					const auto letter = expandedRead[curReadIndex + i];
					insertions[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
					insertions[refGenIndex + i].second.insert(name);
				}
				else if (insertions[refGenIndex + i].second.find(name) == insertions[refGenIndex + i].second.
				                                                                                                end()) {
					//Due to the difference in indexing, we may have the same sequence counted twice for the same
					insertions[refGenIndex + i].first.increase('-');
				}
			}

			if (isInsertion) insertionsFound += iter->second;
			curReadIndex += iter->second;
		}

		//If the section before the end clipping is not an insertion one, the remaining string still needs to be read
		if (readFromIndex != curReadIndex)
			noInsertionsRead += expandedRead.substr(
				readFromIndex, curReadIndex - readFromIndex);

		// The BAM library stores the position for the aligned bases
		startingPos[readSDStart].insert(std::make_pair(noInsertionsRead, name));
	}

	bam_destroy1(b);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
	bam_hdr_destroy(header);
	sam_close(in);

	// while (reader.GetNextAlignment(alignment)) {
	// 	//Check whether the sequence has been aligned to the reference genome
	// 	if (!alignment.IsMapped()) continue;
	//
	// 	if (alignment.Name == "NC-001416-mutated_12202_aligned_166_R_24_12985_15") {
	// 		cout << 1;
	// 	}
	//
	// 	// Clip the cigar in order to remove the clippers to be able to apply the cigar string to the expandedRead directly
	// 	auto cigarExpanded = alignment.CigarData;
	// 	if (cigarExpanded.begin()->Type == 'S' || cigarExpanded.begin()->Type == 'H')
	// 		cigarExpanded.erase(
	// 			cigarExpanded.begin());
	// 	if ((cigarExpanded.end() - 1)->Type == 'S' || (cigarExpanded.end() - 1)->Type == 'H')
	// 		cigarExpanded.erase(
	// 			cigarExpanded.end() - 1);
	//
	// 	//AlignedBases contains the string that is clipped and has '-' inserted at the position of D in the cigar string
	// 	string expandedRead = alignment.AlignedBases;
	// 	//Aligned position for the string with cut out insertions (for the direct substitution and deletion analysis)
	// 	size_t readSDStart = alignment.Position;
	// 	size_t startPos = alignment.Position;
	//
	// 	// If we analyze parts of the string that are not within the window, we can find ourselves in the situation when indices are analyzed more than once with different content
	// 	// check whether the starting index is within the window
	// 	// if so, iterate over the cigar string clipping it in a way that would cover only the part of the read within the window
	// 	// store the index of the last analyzed element in the Cigar string
	// 	// if the starting index is not within the window check the array with stored indices in order to figure out what the position that should be started with and remove the index
	//
	// 	size_t readOffset = 0;
	// 	size_t cigarOffset = 0;
	// 	if (cigarIndices.find(alignment.Name) != cigarIndices.end()) {
	// 		auto aux = cigarIndices[alignment.Name];
	// 		const auto end = next(cigarExpanded.begin(), get<0>(cigarIndices[alignment.Name]));
	// 		cigarExpanded.erase(cigarExpanded.begin(), end);
	// 		cigarExpanded.begin()->Length -= get<1>(cigarIndices[alignment.Name]);
	// 		expandedRead = expandedRead.substr(get<2>(cigarIndices[alignment.Name]));
	// 		cigarOffset = get<0>(cigarIndices[alignment.Name]);
	// 		readOffset = get<2>(cigarIndices[alignment.Name]);
	// 		readSDStart = from;
	// 		startPos = from;
	//
	// 		cigarIndices.erase(alignment.Name);
	// 	}
	// 	else if (startPos < from) continue;
	// 	//Since the insertions are cut out, if the first non-clipped part of the string requires an insertion
	// 	//the starting index should take that into consideration
	// 	if (cigarExpanded.begin()->Type == 'I') readSDStart += cigarExpanded.begin()->Length;
	// 	//Base index for the insertion - should take care of the potential clipping
	//
	// 	size_t insertionsFound = 0;
	// 	string noInsertionsRead;
	// 	size_t readFromIndex = 0;
	// 	size_t curReadIndex = 0;
	//
	// 	for (auto iter = cigarExpanded.begin(); iter != cigarExpanded.end(); iter++) {
	// 		const size_t refGenIndex = startPos - insertionsFound + curReadIndex;
	// 		const bool isInsertion = iter->Type == 'I';
	//
	// 		//if we get out of the boundaries of the defined window
	// 		if (refGenIndex + iter->Length >= to) {
	// 			//We can cover the insertions even if they get out of the window as it does not increase the indices for S and D
	// 			const size_t left = to - refGenIndex;
	// 			const size_t end = isInsertion ? iter->Length : left;
	// 			auto* curMap = &insertions;
	//
	// 			if (isInsertion) {
	// 				noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
	// 				readFromIndex = curReadIndex + iter->Length;
	// 			}
	//
	// 			for (size_t i = 0; i != end; i++) {
	// 				if (i == left) {
	// 					curMap = &insertionsOOB;
	// 				}
	//
	// 				if (isInsertion) {
	// 					(*curMap)[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
	// 					(*curMap)[refGenIndex + i].second.insert(alignment.Name);
	// 				}
	// 				else if ((*curMap)[refGenIndex + i].second.find(alignment.Name) == (*curMap)[refGenIndex + i].
	// 				                                                                   second.
	// 				                                                                   end()) {
	// 					(*curMap)[refGenIndex + i].first.increase('-');
	// 				}
	// 			}
	//
	// 			if (!isInsertion) {
	// 				//We want to store the position in the iteration string, length covered and position in the expanded string
	// 				//in order to start the analysis of the string within the next window quicker
	// 				cigarIndices[alignment.Name] = make_tuple(distance(cigarExpanded.begin(), iter) + cigarOffset, left,
	// 				                                          curReadIndex + left + readOffset);
	// 				curReadIndex += left;
	// 				break;
	// 			}
	//
	// 			insertionsFound += iter->Length;
	// 			curReadIndex += iter->Length;
	// 			continue;
	// 		}
	//
	//
	// 		//TODO ensure that only the parts of the string within the window are analyzed in order to avoid the same indices to be analyzed several times and other errors
	// 		if (isInsertion) {
	// 			noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
	// 			readFromIndex = curReadIndex + iter->Length;
	// 		}
	//
	// 		for (size_t i = 0; i != iter->Length; i++) {
	// 			if (refGenIndex + i == 1056) {
	// 				cout << 2;
	// 			}
	// 			if (isInsertion) {
	// 				const auto letter = expandedRead[curReadIndex + i];
	// 				insertions[refGenIndex + i].first.increase(expandedRead[curReadIndex + i]);
	// 				insertions[refGenIndex + i].second.insert(alignment.Name);
	// 			}
	// 			else if (insertions[refGenIndex + i].second.find(alignment.Name) == insertions[refGenIndex + i].second.
	// 			                                                                                                end()) {
	// 				//Due to the difference in indexing, we may have the same sequence counted twice for the same
	// 				insertions[refGenIndex + i].first.increase('-');
	// 			}
	// 		}
	//
	// 		if (isInsertion) insertionsFound += iter->Length;
	// 		curReadIndex += iter->Length;
	// 	}
	//
	// 	//If the section before the end clipping is not an insertion one, the remaining string still needs to be read
	// 	if (readFromIndex != curReadIndex)
	// 		noInsertionsRead += expandedRead.substr(
	// 			readFromIndex, curReadIndex - readFromIndex);
	//
	// 	// The BAM library stores the position for the aligned bases
	// 	startingPos[readSDStart].insert(std::make_pair(noInsertionsRead, alignment.Name));
	// }
	//
	// set<string> aux;
	// for (const auto& elem : cigarIndices) {
	// 	aux.insert(elem.first);
	// }
	// set<string> readNames = {
	// 	"NC-001416-mutated_17626_aligned_95_R_6_15515_38",
	// 	"NC-001416-mutated_18528_aligned_57_F_13_18082_35",
	// 	"NC-001416-mutated_18774_aligned_77_R_23_14268_3",
	// 	"NC-001416-mutated_19999_aligned_53_F_0_13189_23",
	// 	"NC-001416-mutated_20229_aligned_185_R_20_12433_28",
	// 	"NC-001416-mutated_21158_aligned_63_R_15_11549_42",
	// 	"NC-001416-mutated_21342_aligned_69_F_34_11013_26",
	// 	"NC-001416-mutated_21404_aligned_90_F_25_12182_11",
	// 	"NC-001416-mutated_21446_aligned_156_F_2216_10800_1619",
	// 	"NC-001416-mutated_22101_aligned_16_F_45_11436_23",
	// 	"NC-001416-mutated_22305_aligned_65_F_14_15455_115",
	// 	"NC-001416-mutated_22338_aligned_135_R_18_14296_39",
	// 	"NC-001416-mutated_22914_aligned_76_F_31_12906_21",
	// 	"NC-001416-mutated_22965_aligned_102_R_11_12148_26",
	// 	"NC-001416-mutated_23253_aligned_39_F_34_13492_1",
	// 	"NC-001416-mutated_23432_aligned_182_F_20_9583_19",
	// 	"NC-001416-mutated_24224_aligned_34_F_14_13467_29",
	// 	"NC-001416-mutated_24313_aligned_36_R_30_8412_28",
	// 	"NC-001416-mutated_24467_aligned_218_F_19_8216_26",
	// 	"NC-001416-mutated_24551_aligned_217_F_2_12796_10",
	// 	"NC-001416-mutated_25007_aligned_75_R_13_9294_4",
	// 	"NC-001416-mutated_25327_aligned_161_R_23_14036_23",
	// 	"NC-001416-mutated_25664_aligned_96_R_34_8459_25",
	// 	"NC-001416-mutated_25922_aligned_159_F_17_10884_23",
	// 	"NC-001416-mutated_26192_aligned_119_F_44_13800_20",
	// 	"NC-001416-mutated_26585_aligned_27_R_27_9658_30",
	// 	"NC-001416-mutated_26897_aligned_176_F_0_10728_28",
	// 	"NC-001416-mutated_27077_aligned_107_R_40_6977_30",
	// 	"NC-001416-mutated_27405_aligned_41_R_45_7963_44",
	// 	"NC-001416-mutated_27597_aligned_89_R_32_9958_21",
	// 	"NC-001416-mutated_27633_aligned_197_F_2113_4698_325",
	// 	"NC-001416-mutated_27733_aligned_233_F_88_15361_9",
	// 	"NC-001416-mutated_27790_aligned_64_F_18_4999_19",
	// 	"NC-001416-mutated_28131_aligned_59_F_1593_4753_2108",
	// 	"NC-001416-mutated_28312_aligned_130_F_28_13495_67",
	// 	"NC-001416-mutated_28323_aligned_229_R_28_6853_6",
	// 	"NC-001416-mutated_28445_aligned_216_R_1_10177_12",
	// 	"NC-001416-mutated_28703_aligned_241_R_23_6616_32",
	// 	"NC-001416-mutated_28978_aligned_180_R_32_5854_80",
	// 	"NC-001416-mutated_29034_aligned_238_F_35_5347_24",
	// 	"NC-001416-mutated_29227_aligned_43_F_11_9183_1",
	// 	"NC-001416-mutated_29280_aligned_219_F_19_7706_38",
	// 	"NC-001416-mutated_29344_aligned_78_F_25_15113_1",
	// 	"NC-001416-mutated_29594_aligned_15_F_25_9722_24",
	// 	"NC-001416-mutated_29727_aligned_171_R_45_7125_3",
	// 	"NC-001416-mutated_29800_aligned_170_F_70_3472_33",
	// 	"NC-001416-mutated_29969_aligned_132_F_24_10910_22",
	// 	"NC-001416-mutated_30068_aligned_24_F_1170_2960_84",
	// 	"NC-001416-mutated_30246_aligned_162_F_37_12629_26",
	// 	"NC-001416-mutated_30414_aligned_220_R_46_14056_5",
	// 	"NC-001416-mutated_30419_aligned_189_F_25_7095_24",
	// 	"NC-001416-mutated_30446_aligned_160_R_6_10718_22",
	// 	"NC-001416-mutated_30691_aligned_45_F_38_9635_29",
	// 	"NC-001416-mutated_30710_aligned_193_F_42_8230_25",
	// 	"NC-001416-mutated_30777_aligned_85_F_20_9894_11",
	// 	"NC-001416-mutated_30958_aligned_141_R_63_13311_1",
	// 	"NC-001416-mutated_31504_aligned_133_R_2801_8888_1569",
	// 	"NC-001416-mutated_31668_aligned_237_R_31_12331_39",
	// 	"NC-001416-mutated_31840_aligned_227_F_6_5014_12"
	// };
	// set<string> diff;
	// set_difference(aux.begin(), aux.end(), readNames.begin(), readNames.end(), std::inserter(diff, diff.begin()));
	return {startingPos, insertions, insertionsOOB};
}
