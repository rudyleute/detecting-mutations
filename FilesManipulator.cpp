#include "FilesManipulator.h"

#include <iostream>

#include <filesystem>
#include <fstream>
#include <list>
#include <boost/icl/interval_map.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <map>
#include <set>
#include <utility>
#include <string>

#include "Comparator.h"

using namespace std;
using FM = FilesManipulator;

std::map<string, tuple<size_t, size_t, size_t>> FM::cigarIndices;

string FM::getRefGen(const string& fileName) {
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

string FM::formFullPath(const string& fileName) {
	return filesystem::current_path().parent_path().string() + '/' + fileName;
}

CigarString FM::getCigarString(const bam1_t* b) {
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

MutationsVCF FM::getVCFInsertions(const string& ref, const string& alt, const size_t& pos) {
	MutationsVCF indels;

	char diff = 0;
	for (char c : ref) diff ^= c;
	for (char c : alt) diff ^= c;

	indels[pos] = make_pair(diff, 'I');
	return indels;
}

string FM::getRead(const bam1_t* b) {
	auto* seq = bam_get_seq(b);
	int32_t seqLen = b->core.l_qseq;

	string read(seqLen, 'a');
	for (int i = 0; i < seqLen; i++) {
		switch (bam_seqi(seq, i)) {
		case 1: read[i] = 'A';
			break;
		case 2: read[i] = 'C';
			break;
		case 4: read[i] = 'G';
			break;
		case 8: read[i] = 'T';
			break;
		case 15: read[i] = 'N';
			break;
		default: read[i] = '?';
			break;
		}
	}

	return read;
}

string FM::getExpandedRead(string read, CigarString& cigar) {
	string expandedRead;
	if (cigar.begin()->first == 'S' || cigar.begin()->first == 'H') {
		read = read.substr(cigar.begin()->second);
		cigar.erase(cigar.begin());
	}

	auto cigarLast = next(cigar.end(), -1);
	if (cigarLast->first == 'S' || cigarLast->first == 'H') {
		read = read.substr(0, read.size() - cigarLast->second);
		cigar.erase(cigarLast);
	}

	size_t readFrom = 0;
	size_t curInd = 0;
	for (const auto& cigarElem : cigar) {
		if (cigarElem.first == 'D') {
			expandedRead += read.substr(readFrom, curInd - readFrom);
			readFrom = curInd;

			expandedRead += string(cigarElem.second, '-');
		} else curInd += cigarElem.second;
	}

	if (curInd != readFrom) expandedRead += read.substr(readFrom, curInd - readFrom);

	return expandedRead;
}

void FM::saveToCsv(const string& geneName, CompRes& errors) {
	std::map<size_t, string> indices;

	for (const auto &aux: errors.diffInVCF) {
		ostringstream str;
		str << "Missed, " << std::get<0>(aux) << ", " << std::get<2>(aux) << ", " << std::get<1>(aux);
		for (const auto aux2: std::get<3>(aux).getCounters()) str << ", " << aux2;
		indices[std::get<0>(aux)] = str.str();
	}

	for (const auto &aux: errors.diffInCust) {
		ostringstream str;
		str << "Additional, " << std::get<0>(aux) << ", " << std::get<2>(aux) << ", " << std::get<1>(aux);
		for (const auto aux2: std::get<3>(aux).getCounters()) str << ", " << aux2;
		indices[std::get<0>(aux)] = str.str();
	}

	for (const auto &aux: errors.errors) {
		ostringstream str;
		str << "Error, " << std::get<0>(aux) << ", " << std::get<2>(aux) << ", " << std::get<1>(aux);
		for (const auto aux2: std::get<5>(aux).getCounters()) str << ", " << aux2;
		str << ", " << std::get<3>(aux) << ", " << std::get<4>(aux);
		indices[std::get<0>(aux)] = str.str();
	}

	ofstream csvOut(FM::formFullPath(geneName + ".csv"));

	csvOut << "Type, Index, Action, Symbol, ";
	for (const auto & [fst, snd]: nucleoMapping) csvOut << fst << ", ";
	csvOut << "Expected Nucleo, Expected Action" << endl;

	for (const auto & [fst, snd]: indices) csvOut << snd << endl;

	csvOut.close();
}

MutationsVCF FM::readFreeBayesVCF(const string& fileName) {
	htsFile* fp = bcf_open(fileName.c_str(), "r");
	if (!fp) {
		cerr << "Failed to open the file " << fileName << endl;
		throw runtime_error("Failed to open the file " + fileName);
	}

	bcf_hdr_t* hdr = bcf_hdr_read(fp);
	bcf1_t* rec = bcf_init();

	MutationsVCF mutations;
	while (bcf_read(fp, hdr, rec) == 0) {
		bcf_unpack(rec, BCF_UN_STR);

		int pos = rec->pos;

		// REF allele is first allele
		const string ref = rec->d.allele[0];
		if (rec->n_allele < 2) {
			std::cerr << "No ALT alleles at pos " << pos << "\n";
			continue;
		}

		const string alt = rec->d.allele[1];
		size_t refLen = ref.size();
		size_t altLen = alt.size();

		std::string variantType;
		if (refLen == 1 && altLen == 1) mutations[pos] = std::make_pair(alt[0], 'X');
		else {
			size_t lenDiff = (refLen > altLen) ? (refLen - altLen) : (altLen - refLen);

			if (lenDiff >= 2) mutations[pos + 1] = std::make_pair('U', 'C');
			else if (refLen > altLen) mutations[pos + 1] = std::make_pair('-', 'D');
			else if (refLen < altLen) mutations.merge(getVCFInsertions(ref, alt, pos + 1));
		}
	}

	// Cleanup
	bcf_destroy(rec);
	bcf_hdr_destroy(hdr);
	bcf_close(fp);

	return mutations;
}

size_t FM::getRefGenLength(const string& fileName) {
	samFile* in = sam_open(fileName.c_str(), "r");
	if (!in) {
		cerr << "Error opening file " << fileName << endl;
		throw runtime_error("Error opening file " + fileName);
	}

	bam_hdr_t* header = sam_hdr_read(in);
	const size_t refLen = header->target_len[0];

	bam_hdr_destroy(header);
	sam_close(in);

	return refLen;
}

string FM::getRefGenName(const string& fileName) {
	samFile* in = sam_open(fileName.c_str(), "r");
	if (!in) {
		cerr << "Error opening file " << fileName << endl;
		throw runtime_error("Error opening file " + fileName);
	}

	bam_hdr_t* header = sam_hdr_read(in);
	const string refName = header->target_name[0];

	bam_hdr_destroy(header);
	sam_close(in);

	return refName;
}

AlignmentMaps FM::getAlignments(
	const string& fileName,
	const size_t& from,
	const size_t& to,
	const string& refName,
	const InsertionMap& prevIterInsertions
) {
	samFile* in = sam_open(fileName.c_str(), "r");
	bam_hdr_t* header = sam_hdr_read(in);
	hts_idx_t* idx = sam_index_load(in, fileName.c_str());
	hts_itr_t* iter = sam_itr_querys(idx, header,
	                                 (refName + ":" + std::to_string(from) + "-" + std::to_string(to)).c_str());
	bam1_t* b = bam_init1();

	std::map<size_t, set<pair<string, string>>> startingPos;
	Insertions insertions = Insertions(prevIterInsertions);

	while (sam_itr_next(in, iter, b) >= 0) {
		//Check whether the sequence has been aligned to the reference genome
		if (b->core.flag & BAM_FUNMAP) continue;

		size_t pos = b->core.pos;
		string name = bam_get_qname(b);
		CigarString cigarExpanded = getCigarString(b);

		const string read = getRead(b);
		string expandedRead = getExpandedRead(read, cigarExpanded);

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
			const auto end = next(cigarExpanded.begin(), get<0>(cigarIndices[name]));
			cigarExpanded.erase(cigarExpanded.begin(), end);
			cigarExpanded.begin()->second -= get<1>(cigarIndices[name]);
			expandedRead = expandedRead.substr(get<2>(cigarIndices[name]));
			cigarOffset = get<0>(cigarIndices[name]);
			readOffset = get<2>(cigarIndices[name]);
			readSDStart = from;
			startPos = from;

			cigarIndices.erase(name);
		} else if (startPos < from) continue;
		//Since the insertions are cut out, if the first non-clipped part of the string requires an insertion
		//the starting index should take that into consideration
		if (cigarExpanded.begin()->first == 'I') readSDStart += cigarExpanded.begin()->second;
		//Base index for the insertion - should take care of the potential clipping

		size_t insertionsFound = 0;
		string noInsertionsRead;
		size_t readFromIndex = 0;
		size_t curReadIndex = 0;

		insertions.setRead(expandedRead, name);
		for (auto cigarIter = cigarExpanded.begin(); cigarIter != cigarExpanded.end(); cigarIter++) {
			const size_t refGenIndex = startPos - insertionsFound + curReadIndex;
			const bool isInsertion = cigarIter->first == 'I';

			//if we get out of the boundaries of the defined window
			if (refGenIndex + cigarIter->second >= to) {
				//We can cover the insertions even if they get out of the window as it does not increase the indices for S and D
				const size_t left = to - refGenIndex;
				const size_t end = isInsertion ? cigarIter->first : left;

				if (isInsertion) {
					noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
					readFromIndex = curReadIndex + cigarIter->second;
				}

				insertions.addInsertion(refGenIndex, curReadIndex, end, left, isInsertion);
				if (!isInsertion) {
					//We want to store the position in the iteration string, length covered and position in the expanded string
					//in order to start the analysis of the string within the next window quicker
					cigarIndices[name] = make_tuple(distance(cigarExpanded.begin(), cigarIter) + cigarOffset, left,
					                                curReadIndex + left + readOffset);
					curReadIndex += left;
					break;
				}

				insertionsFound += cigarIter->second;
				curReadIndex += cigarIter->second;
				continue;
			}

			if (isInsertion) {
				noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);
				readFromIndex = curReadIndex + cigarIter->second;
			}
			insertions.addInsertion(refGenIndex, curReadIndex, cigarIter->second, cigarIter->second, isInsertion);

			if (isInsertion) insertionsFound += cigarIter->second;
			curReadIndex += cigarIter->second;
		}

		//If the section before the end clipping is not an insertion one, the remaining string still needs to be read
		if (readFromIndex != curReadIndex)
			noInsertionsRead += expandedRead.substr(readFromIndex, curReadIndex - readFromIndex);

		// The BAM library stores the position for the aligned bases
		startingPos[readSDStart].insert(make_pair(noInsertionsRead, name));
	}

	bam_destroy1(b);
	hts_itr_destroy(iter);
	hts_idx_destroy(idx);
	bam_hdr_destroy(header);
	sam_close(in);

	return {startingPos, insertions};
}
