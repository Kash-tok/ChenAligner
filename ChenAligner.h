#pragma once

#include <string>
#include <chrono>
#include <stdexcept>
#include <vector>
#include <map>
#include <queue>
#include <list>
#include <stack>
#include <tuple>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>
#include <mutex>
#include <condition_variable>
#include <cmath>
#include <atomic>
#include <algorithm>

/************************************************************
 ** Part Segment Couple *************************************
 * This part deals with non-overlapping segments ************
 * It is used for 'coverage length' and 'watermark ranges' **
 ***********************************************************/

struct SCNode
{
	unsigned int l, r;
	SCNode* next = nullptr;

	SCNode(unsigned int lvalue, unsigned int rvalue) : l(lvalue), r(rvalue) {}
};

class SegmentCouple
{
private:
	SCNode* m_begin = nullptr;
	unsigned int m_coverlength = 0;
	unsigned int minl, maxr;

public:
	SegmentCouple() = default;
	~SegmentCouple() = default;

	void insert_noOverlap(unsigned int segl, unsigned int segr);	// keep all segments that don't overlap
	void insert(unsigned int segl, unsigned int segr);				// combine segments that overlap
	void print_noOverlap();
	void print();
	unsigned int coverlength() const { return m_coverlength; }
	std::vector<unsigned int> best_joined();
	void best_joined(std::ostream& os);
	std::vector<std::pair<unsigned int, unsigned int>> get_all_segments() const;
	std::vector<bool> get_individual_cover(unsigned int l, unsigned int r) const;
};

class SegmentCouples
{
private:
	std::map<unsigned int, SegmentCouple> storage;
	unsigned int max_coverlength = 0;
	std::vector<std::pair<unsigned int, SegmentCouple>> good_matches;

	friend class ChenAlignerCore;

public:
	void insert(unsigned int beginid, unsigned int endid, unsigned int segl, unsigned int segr);
	void insert(std::vector<unsigned int> id, unsigned int segl, unsigned int segr);
	void insert(const std::vector<std::pair<unsigned int, unsigned int>>& cp, unsigned int segl, unsigned int segr);
	unsigned int maxcoverlength() { return max_coverlength; }
	std::vector<std::pair<unsigned int, SegmentCouple>> show_matches(unsigned int density_diff, std::vector<unsigned int> id_convertor = {});
};


/******************************************
 ** Part 'Judger' *************************
 *****************************************/

#define USE_STAT_DATA
const int STAT_DATA[12][21] = { {15,21,77,93,94,92,92,94,76,48,0,0,0,0,0,0,0,0,0,0,0},{14,21,79,94,93,95,93,91,91,92,84,85,53,44,0,0,0,0,0,0,0},{14,22,92,95,93,95,96,94,93,88,85,85,88,88,81,84,87,92,0,0,0},{13,24,92,94,95,97,94,93,93,97,93,92,97,76,79,63,43,44,0,0,0},{13,25,91,94,96,97,93,99,93,94,94,92,97,93,92,89,88,93,74,79,76},{13,28,88,94,95,97,98,96,94,95,96,93,94,91,92,91,92,91,90,87,77},{13,30,78,94,97,98,98,98,96,102,96,94,101,95,95,95,92,94,95,96,92},{13,34,78,94,96,101,101,97,97,96,96,97,95,96,95,97,94,96,98,100,101},{12,38,47,96,96,100,97,99,100,100,97,98,98,99,100,100,101,99,101,100,99},{12,43,54,94,96,99,101,98,102,99,99,102,102,100,105,101,101,101,100,98,100},{12,49,59,95,101,102,103,104,105,104,104,105,104,106,104,107,103,102,102,99,91},{13,56,67,96,99,103,105,107,111,110,112,107,113,109,111,105,110,92,93,93,0} };

/* IN THIS PART, ALL THE 'LOG' ARE ACTURALLY '-LOG' */

/* m out of n */
inline double logC(unsigned int m, unsigned int n)
{
	return  std::lgamma(m + 1) + std::lgamma(n - m + 1) - std::lgamma(n + 1);
}

#ifndef USE_STAT_DATA
/* This class calculates the probability of having a random string that has <=err hamming distance with a given string of length 'len' */
class ErrorRateConventor
{
private:
	static std::atomic_bool m_loglikInited;
	static std::mutex m_Mutex;
	static std::vector<double> m_loglik;
	unsigned int m_alphabetSize;
	unsigned int m_maxLen;
	unsigned int m_standardError;

	double stLogP;

public:
	ErrorRateConventor(unsigned int maxlen, unsigned int alphabet_size, unsigned int standard_error = 0);

	double logP(const unsigned int& len, const unsigned int& err) const;
	bool goodMatch(const unsigned int& len, const unsigned int& err) const;
};
#else
/* This class calculates the probability of having a random string that has <=err hamming distance with a given string of length 'len' */
class ErrorRateConventor
{
private:
	static std::atomic_bool m_loglikInited;
	static std::mutex m_Mutex;
	static std::vector<double> m_minl;

public:
	ErrorRateConventor(unsigned int maxlen, unsigned int alphabet_size, unsigned int standard_error = 0);

	bool goodMatch(const unsigned int& len, const unsigned int& err) const;
};
#endif // !USE_STAT_DATA


/************************************
 ** Core of the aligner *************
 * For developers *******************
 ***********************************/

constexpr bool SHOW_OUTPUT_enlengthenSearch_ignoreErasure = false;
constexpr bool SHOW_OUTPUT_print_segment = false;
constexpr bool SHOW_OUTPUT_enlengthen_procedure = false;
// #define SHOW_OUTPUT_substring_history
constexpr int INIT_TYPE = 0;
constexpr unsigned int confidence_requirement = 1;


/* The substring that has limited Hamming distance with any substring of watermarks.
* See report for details.
*/
struct MatchedSubstring
{
	unsigned int node_num;		// the node at which link must be used
	unsigned int pattern_pos;	// the position of the end+1 of the substring in pattern string
	unsigned int matched_len;	// the length of the substring
	unsigned int error_number;

	bool operator < (const MatchedSubstring& other) const
	{
		return matched_len > other.matched_len;
	}
};

struct GSANode
{
	// standard part
	int len = 0;	// length of this substring = depth on the trie
	int link = 0;	// index of the longest suffix that appealled more times than this string
	std::map<char, unsigned int> next;

	// extra part
	unsigned int pos = 0;	// the end position +1 of this string when it _first_ appeals in the main string (empty string.pos = 0)
	unsigned int host = 0;	// the index of the _first_ main string that contains this substring
	std::vector<unsigned int> inverse_link;		// inverse of link
	//UNUSED int parent = 0;// the _least_ deep node on the tree that has 'parallel' link with this node 
							// or the end+1 of the prefix that is discarded after using link
							// if linked to 0, the parent will be itself

	// for cache storage
	void print(std::ostream& os) const;
	void read(std::istream& is);
};

enum class OUTFILE_FORMAT
{
	BEST, ALL
};

/*
* This class uses the traditional methods in segment tree, but don't have lazy tags.
* In update, a binary-tree-like structure that maintains the minimum number of its children is maintained.
* In query, the maximum number through the searching pathway is returned.
* The structure guarantees that the pathway is unimodal.
* Overall, the sturcture maintains the maximum number that one placed had during the update.
*/
class SegmentTree
{
private:
	bool m_updateHappened = true;
	std::vector<unsigned int> m_STMin;
	constexpr unsigned int stlchild(unsigned int x) { return x << 1; }
	constexpr unsigned int strchild(unsigned int x) { return (x << 1) + 1; }
	constexpr unsigned int stparent(unsigned int x) { return x >> 1; }

public:
	void Update(const unsigned int& left, const unsigned int& right, const unsigned int& value, const unsigned int& segleft, const unsigned int& segright, const unsigned int& node);
	unsigned int Query(const unsigned int& x, const unsigned int& segleft, const unsigned int& segright, const unsigned int& node);
	SegmentTree(size_t size);
	void updateTagReset() { m_updateHappened = false; }
	bool updateHappened() const { return m_updateHappened; }
};

std::string read2filename(const std::string& dna_sequence);

class ChenAlignerCore
{
private:
	std::string channelTypeName = "ERASE0.050000";

	/// sotrage
	bool m_IsEraseChannel;
	char m_EraseMarker;
	unsigned int m_TerminationErrorCount;
	unsigned int m_CharacterNum;	// size of dictionary
	unsigned int m_StringCount = 0;	// total number of host strings
	std::set<char> m_Alphabet;
	std::vector<GSANode> m_Nodes;
	std::vector<std::pair<unsigned int, unsigned int>> m_HostRangeTrie;		// built from trie, for each node, store <min, max> of the hosts sorted dictionarily. Don't store to cache.
	std::vector<SegmentCouple> m_HostRange;		// Don't store to cache (Host = Watermark in the report)
	std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_HostRangeCache;
	std::vector<unsigned int> m_HostMapping;	// map the order to id of host string
	std::vector<unsigned int> m_NodeBefore;		// inverse of next
	std::vector<int> m_DeepestLink;				// when there are chain of erase markers, the deepest link of every node on the chain
	std::vector<unsigned int> m_EraseChainLen;	// the number of continus erasure markers until this one

	/// working variables
	unsigned int m_LastNode = 0;	// index of the last node visited or inserted
	unsigned int m_StringId;		// identity of the main string being added to trie
	bool m_HaveBuilt = false;		// adding new strings is only legal before building the SAM, when this value is still "false"
	unsigned int m_UnclonedLength;

	/// Methods
	/// Trie
	inline void insertTrie(const char& c);			// append one char to the end of the string stored in trie
	void calHostRange(const unsigned int& node);	// calculate the range of watermarks each node belong to. It is not usually useful since watermarks are not similar
	// it is still useful when erase/error rate is high, and the matched segments are short

	/// SAM (these functions are standard way of building a suffix automation
	inline unsigned int insertSAM(const unsigned int& last, const char& c);
	inline bool isNotCloned(const unsigned int& node) const;	// if index larger than tire size, then it's cloned. cloned string shouldn't be counted
	// similar to exact matching in algorithm, but to be used in full match
	inline std::vector<unsigned int> substringExtension(const unsigned int& start_node, const unsigned int& start_len, const unsigned int& max_successor,
		const std::string& str, const unsigned int start_pos, unsigned int& best_len, unsigned int& best_node, unsigned int& best_pos) const;
	inline std::vector<unsigned int> substringExtensionWithErasure(const unsigned int& start_node, const unsigned int& start_len, const unsigned int& max_successor,
		const std::string& str, const unsigned int start_pos, unsigned int& best_len, unsigned int& best_node, unsigned int& best_pos, const char& ignored_char) const;

	void buildHostRangeTire();	// To be called before building SAM, host range in trie
	void buildSAM();
	void buildInverseLink();	// After building SAM. The inversed links should form a tree.
	void buildInverseNext();	// After building SAM. Calculate the direct 'parent' of each node. Exactly one layer/len above.
	void calFullHostRange(const unsigned int& node);	// inner function of buildHostRangeFull
	void buildHostRangeFull();	// After building inverse link. Use segment couples to cover all the hosts.
	void buildDeepestLink();	// After building SAM. Only use if m_IsErasureChannel is true.
	// It also calculates erase chain length

	/// Parallel
	std::vector<std::vector<unsigned int>> m_Results;  // the results
	std::vector<std::thread> m_Matchers;	// workers in thread pool
	std::vector<std::string> m_TaskReads;	// reads
	unsigned int m_TasksTaken;		// workers pick tasks before they do
	std::mutex m_Mutex;						// traditional locking system
	std::condition_variable m_ConVar;		// traditional locking system
	unsigned int m_sub_err;
	unsigned int m_all_err;
	unsigned int m_err_dif;

	void workerThread();	// a worker, as its name says
	void initThreadPool(unsigned int thread_num = 10);	// create a thread pool. Automatically called by parallelMatch if not there.

	inline void print_host_pairs(const std::vector<std::pair<unsigned int, unsigned int>> data, std::ostream& os) const;

	/// erase channel special
	std::vector<std::pair<unsigned int, SegmentCouple>> alignReadWithErasure(const std::string& str, const char& erased_char, const unsigned int& nature_sub_err, const unsigned int& nature_other_err, const unsigned int& err_dif) const;

public:
	ChenAlignerCore(unsigned int dict_size);
	ChenAlignerCore(unsigned int dict_size, char erase_marker, unsigned int termination_err_cnt = 5);
	~ChenAlignerCore() = default;

	/// prepare (call these functions in order)
	bool addString(const std::string& str, const unsigned int& id);		// add a new string to trie; if already built, return false; id >= 0
	bool addString(const char* str, const unsigned int& id);
	void build();

	/// use

	/* Exact Match Ot(n) Om(1) (usually not used) *
	* input: the string 'str' to be matched
	* output: <id, pos> where 'pos' is the place where 'str' starts in reference string of identity 'id'
	* failure: return <0,0>
	*/
	std::pair<unsigned int, unsigned int> oldFunction_exactAlignment(const std::string& str);

	/* Exact Common Substring Ot(n) Om(1) (usually not used) *
	* input: the string 'str'
	* output: <id, pos, st, len>
	* 	pos: start position of substring in reference string
	*	st: start position in 'str'
	*	len: matched length
	* substring: the longest one
	*/
	std::tuple<unsigned int, unsigned int, unsigned int, unsigned int> oldFunction_exactCommonSubstring(const std::string& str);

	/* Substitution Match Ot(n) Om(E*n^2) (usually not used) *
	* input: the string 'str', maximum allowed mismatch, output vector length
	* output: vector<id, start_pos_host, start_pos_pat, len, score(mismatchNo.)> sorted by score/len
	*
	* The beginning and end of the result substring is always a match
	*/
	std::vector< std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int>> oldFunction_alignWithSubstitution(const std::string& str, const unsigned int& maxE);

	/* Full Match (the core of ChenAligner)
	* by all means, this aligner was initially designed for flipping channels, rather than erasing channels
	* input:	read 'str',
	*			equivalent maximum substitution count in a string of length of 'str' (for substring matching),
	*			expected deletion density as equivalent substring length in a string of length of 'str' (for full string matching),
	*			maximum absolute density difference compared to the best fit as equivalent length (for output size)
	* output:	std::vector<std::pair<host_index, SegmentCouple>>
	*/
	std::vector<std::pair<unsigned int, SegmentCouple>> alignRead(const std::string& str, const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif) const;


	/// Parallel
	void parallel_Init(const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif);
	void parallel_addRead(const std::string& str);
	void parallel_addRead(const std::string& str, const unsigned int& sub_err, const unsigned int& all_err, const unsigned int& err_dif);
	void parallel_terminateReadTransfer(const unsigned int& thread_num = 10);  // call this function when all reads has been uploaded
	void parallel_getResults(const std::string& outfile_name, OUTFILE_FORMAT format);
	void parallel_getResults(std::vector<unsigned int>& best_result);  // return the index of the watermark that matches best
	void parallel_getResults(std::vector<std::vector<unsigned int>>& all_results);  // return also the indexes of watermarks that are similarly good


	/// for cache storage
	void cache_storeStructure(std::ostream& os) const;
	void cache_acquireStructure(std::istream& is);
};


// debug functions
void print_segment(const std::string& str, const std::vector<std::pair<unsigned int, unsigned int>>& host_range, const unsigned int& l, const unsigned int& r);


/****************************************
 ** Wrapped Chen's Aligner **************
 * For General Users ********************
 ***************************************/

enum class ReadType
{
	FASTQ,
	FASTA,
	PURE
};

class ChenAligner
{
private:
	static constexpr double INDEL_RATE = 0.4;	// the non-substitution-only part length
	static constexpr double LEN_DIFF_ACCEPTED = 0.1;	// the criteria of second-best match
	static constexpr char LOG_FILE_NAME[] = "log.txt";

	ChenAlignerCore m_Core;

	bool init_finished = false;

	void MatchReadsNormal(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate);
	void MatchReadsParallel(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate, unsigned int thread_count);


public:
	/* Define a ChenAligner object
	* Parameters:
	*   alphabet_size: number of types of characters contained in watermarks. It is 4 by default, representing the {A,G,C,T} basepairs.
	*/
	ChenAligner(unsigned int alphabet_size);
	ChenAligner(unsigned int alphabet_size, char erase_marker);

	/* Get watermarks of text form from a file. MatchReads is not available until this function is called.
	* Parameters:
	*   file_name: name of the file that stores all the watermarks in text (not binary) form (please include suffix like '.txt')
	*/
	void GetWatermarksFromFile(const std::string& file_name);

	/* Get watermarks of cache, which is faster. The cache can be built using DumpWatermark method.
	* Parameters:
	*   file_name: name of the file that stores the cache data. This file is created by dump method. (please include suffix like '.cache')
	*/
	void GetWatermarksFromCache(const std::string& file_name);

	/* Store the built watermark structure in cache, and next time GetWatermarksFromCache method can be used.
	* Parameters:
	*   file_name: name of the binary file to store the data (please include suffix like '.cache')
	*/
	void DumpWatermarkStructureToCache(const std::string& file_name);

	/* Match the reads in a file to watermarks.
	A list of indexes of the id of the watermark that each read belongs to will be calculated.
	If a read can't be matched to a watermark, the corresponding index will be 0.
	* Parameters:
	*   read_type: the format of the file that stores reads. Support .fastq .fasta and pure text files
	*   file_name: name of the file that stores reads in read_type format (please include suffix like '.fastq')
	*   result: this vector will be modified to store the matching results
	*   sub_error_rate: a number between 0 and 1 indicating the desity of information artifically added to each watermark. Larger numbers will slow down the matcher, while smaller numbers will result in lower success rate.
	*   use_parallel: set true to allow the program to use multiple CPU cores for faster matching
	*/
	void MatchReads(ReadType read_type, std::string file_name, std::vector<unsigned int>& result, double sub_error_rate, bool use_parallel = false, unsigned int thread_count = 10);
	void MatchReads(ReadType read_type, std::string file_name, std::string result_filename, double sub_error_rate, unsigned int thread_count = 10);
};