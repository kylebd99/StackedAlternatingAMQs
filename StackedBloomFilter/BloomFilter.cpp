#include "BloomFilter.h"



BloomFilter::BloomFilter(size_t filter_size, int num_hashes, int seed)
{
	seed_ = seed;
	num_hashes_ = num_hashes;
	if (filter_size <= 2) filter_size = 2;
	filter_size_ = filter_size;
	filter_.resize(filter_size, false);
}

BloomFilter::~BloomFilter()
{
}

bool BloomFilter::testElement(long element) {
	num_checks_++;
	size_t hash1 = getHash1(element);
	size_t hash2 = getHash2(element);
	for (unsigned int i = 0; i < num_hashes_; i++) {
		if (filter_[getNthHash(hash1, hash2, i)] == false) {
			return false;
		}
	}
	return true;
}

void BloomFilter::addElement(long element) {
	num_elements_++;
	size_t hash1 = getHash1(element);
	size_t hash2 = getHash2(element);
	for (unsigned int i = 0; i < num_hashes_; i++) {
		filter_[getNthHash(hash1, hash2, i)] = true;
	}
}

size_t BloomFilter::getHash1(long x) {
	return CityHash64WithSeed((char*)&x, 4, 1234567+seed_);
}

size_t BloomFilter::getHash2(long x) {
	return CityHash64WithSeed((char*)&x, 4, 7654321+seed_);
}

// Double hashing strategy recommended in Mitzenmacher Paper.
size_t BloomFilter::getNthHash(size_t hash1, size_t hash2, size_t hash_num) {
	return (hash1 + hash_num * hash2) % filter_size_;
}

size_t BloomFilter::GetSize() {
	return filter_size_;
}

double BloomFilter::GetLoadFactor() {
	size_t load = 0;
	for (size_t i = 0; i < filter_size_;  i++) if (filter_[i] == true) load++;
	return (double)load / (double)filter_size_;
}

size_t BloomFilter::GetNumElements() {
	return num_elements_;
}