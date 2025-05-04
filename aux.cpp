#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <cstdint>
#include <algorithm>

#include <api/BamReader.h>     // BamTools
#include <api/BamAlignment.h>

using namespace std;
using namespace BamTools;

// small helpers
inline int     baseIndex(char b) {
    switch (b) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default : return 4;  // use index 4 for deletions/gaps
    }
}

// Walk the CIGAR ops in `cigar` to find what this read contributes at reference position `pos`.
// Returns (‘-’ for a deletion, or the read base for a match/mismatch).
char findBaseAt(const string &seq,
                const vector<pair<char,uint32_t>> &cigar,
                uint32_t refStart,   // 1-based start of this read on ref
                uint32_t pos)        // 1-based ref position to query
{
    uint32_t rpos = refStart;
    uint32_t qpos = 0;  // index into seq

    for (auto &op : cigar) {
        char type = op.first;
        uint32_t len = op.second;

        if (type == 'M' || type == '=' || type == 'X') {
            if (pos < rpos + len) {
                // it falls within this block
                return seq[qpos + (pos - rpos)];
            }
            rpos += len;
            qpos += len;
        }
        else if (type == 'D') {
            if (pos < rpos + len) {
                // deletion at this ref pos
                return '-';
            }
            rpos += len;
        }
        else if (type == 'I') {
            // insertion: consumes read but not ref
            qpos += len;
        }
        else if (type == 'S' || type == 'H') {
            // clipping: consumes read only
            if (type == 'S') qpos += len;
        }
        else {
            // other CIGAR ops (e.g. N, P) can be handled here if needed
        }
    }

    // if we fall off, treat as no coverage
    return '-';
}

struct Sub { uint32_t pos; char refB, altB; uint32_t altCt, depth; };
struct Del { uint32_t pos; uint32_t delCt, depth; };
struct Ins { uint32_t pos; string seq; };

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: call_variants <in.bam> <ref.fa>\n";
        return 1;
    }
    string bamFile = argv[1];
    string refFasta = argv[2];

    // 1) Load reference sequence (naïve)
    vector<char> reference;
    {
        ifstream fin(refFasta);
        string line;
        while (getline(fin, line)) {
            if (line.empty() || line[0]=='>') continue;
            for (char c : line) {
                reference.push_back(c);
            }
        }
    }
    uint32_t refLen = reference.size();   // 0-based; we'll index 1..refLen

    // 2) Prepare event tables
    vector<vector<uint32_t>> starts(refLen+2), ends(refLen+2);
    vector<vector<string>>   inserts(refLen+2);

    // 3) Read BAM, parse alignments
    BamReader reader;
    if (!reader.Open(bamFile)) {
        cerr << "Error opening BAM " << bamFile << "\n";
        return 1;
    }

    BamAlignment aln;
    struct ReadRec {
        uint32_t                   refStart, refEnd;
        string                     seq;
        vector<pair<char,uint32_t>> cigar;
    };
    vector<ReadRec> reads;
    reads.reserve(1e6);

    while (reader.GetNextAlignment(aln)) {
        if (!aln.IsMapped()) continue;

        // parse CIGAR
        vector<pair<char,uint32_t>> cig;
        for (auto &c : aln.CigarData) {
            cig.emplace_back(c.Type, c.Length);
        }

        // get sequence, RC if needed
        string seq = aln.QueryBases;
        if (aln.IsReverseStrand()) {
            reverse(seq.begin(), seq.end());
            for (auto &c : seq) {
                switch(c) {
                    case 'A': c='T'; break;
                    case 'T': c='A'; break;
                    case 'C': c='G'; break;
                    case 'G': c='C'; break;
                    default: break;
                }
            }
        }

        uint32_t a = aln.Position + 1;           // BamTools gives 0-based internally
        uint32_t b = a + aln.AlignedBases.size() - 1;

        uint32_t rid = reads.size();
        reads.push_back({a, b, seq, cig});
        starts[a].push_back(rid);
        ends[b+1].push_back(rid);

        // record insertion sequences
        uint32_t curRef = a;
        uint32_t curSeq = 0;
        for (auto &op : cig) {
            char t = op.first;
            uint32_t L = op.second;
            if (t=='I') {
                inserts[curRef].push_back(seq.substr(curSeq, L));
                curSeq += L;
            }
            else if (t=='M' || t=='=' || t=='X') {
                curRef += L; curSeq += L;
            }
            else if (t=='D') {
                curRef += L;
            }
            else if (t=='S' || t=='H') {
                if (t=='S') curSeq += L;
            }
        }
    }
    reader.Close();

    // 4) Sweep & call
    array<uint32_t,5> base_counts = {0,0,0,0,0};
    vector<Sub> substitutions;
    vector<Del> deletions;
    vector<Ins> insertions;

    for (uint32_t pos = 1; pos <= refLen; ++pos) {
        // remove ended reads
        for (auto rid : ends[pos]) {
            char b = findBaseAt(reads[rid].seq, reads[rid].cigar,
                                reads[rid].refStart, pos-1);
            base_counts[ baseIndex(b) ]--;
        }
        // add starting reads
        for (auto rid : starts[pos]) {
            char b = findBaseAt(reads[rid].seq, reads[rid].cigar,
                                reads[rid].refStart, pos);
            base_counts[ baseIndex(b) ]++;
        }

        // record insertions at this site
        for (auto &seq : inserts[pos]) {
            if (!seq.empty())
                insertions.push_back({pos, seq});
        }

        // depth filter (excluding deletions)
        uint32_t depth = base_counts[0]
                       + base_counts[1]
                       + base_counts[2]
                       + base_counts[3];
        if (depth < 5) continue;

        // substitutions
        char refb = reference[pos-1];
        uint32_t refCt = base_counts[ baseIndex(refb) ];
        if (depth - refCt > refCt) {
            char bestAlt = refb;
            uint32_t bestCt = 0;
            for (char b : {'A','C','G','T'}) {
                auto ct = base_counts[ baseIndex(b) ];
                if (b!=refb && ct>bestCt) {
                    bestAlt = b; bestCt = ct;
                }
            }
            substitutions.push_back({pos, refb, bestAlt, bestCt, depth});
        }
        // deletions
        uint32_t delCt = base_counts[4];
        if (delCt > 0) {
            deletions.push_back({pos, delCt, depth});
        }
    }

    // 5) Output (tab-delimited)
    cout << "#POS\tREF\tALT\tALT_CT\tDEPTH\n";
    for (auto &s : substitutions) {
        cout << s.pos << "\t"
             << s.refB << "\t"
             << s.altB << "\t"
             << s.altCt << "\t"
             << s.depth << "\n";
    }
    cout << "#DELETIONS: POS\tDEL_CT\tDEPTH\n";
    for (auto &d : deletions) {
        cout << d.pos << "\t"
             << d.delCt << "\t"
             << d.depth << "\n";
    }
    cout << "#INSERTIONS: POS\tSEQ\n";
    for (auto &i : insertions) {
        cout << i.pos << "\t" << i.seq << "\n";
    }

    return 0;
}
