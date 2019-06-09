#pragma once
#include "Common.h"

class BloomFilter
{
private:
	std::vector<bool> filter_;
	size_t getHash1(long element);
	size_t getHash2(long element);
	size_t getNthHash(size_t hash1, size_t hash2, size_t hash_num);

public:
	unsigned int num_hashes_;
	size_t filter_size_;
	size_t num_checks_ = 0;
	size_t num_elements_ = 0;
	int seed_;
	BloomFilter(size_t filter_size, int num_hashes, int seed);
	~BloomFilter();
	bool testElement(long element);
	void addElement(long element);
	size_t GetSize();
	double GetLoadFactor();
	size_t GetNumElements();
};

