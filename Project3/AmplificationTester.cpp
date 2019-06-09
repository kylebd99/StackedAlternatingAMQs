#include "../StackedBloomFilter/BloomFilter.h"
#include "../StackedBloomFilter/CityHash.h"
#include <fstream>
#include <stdlib.h>
#include <random>

// Some thoughts... simple amplification doesn't seem to give a lot of boost in performance. 
// It's not nothing, about a 10-15% savings in FPR, but not exactly groundbreaking. Now, if there were a family in which
// previously successful members of the family could point you towards more successful ones, then this might be more
// tractable. Perhaps look into locality perserving hash functions or just the broader family of hash functions? 
// Further, the stacked filter design might see a more fruitful application of this because it likely has a higher
// variance in filter performance in the first place. Could a "worse" hash function benefit enough from this method
// to make up for the clustering of the negative hashes that would be induced as well?



std::vector<long> generate_ints(int num_elements) {
	std::random_device rd;
	/* Random number generator */
	std::default_random_engine generator(rd());
	/* Distribution on which to apply the generator */
	std::uniform_int_distribution<long> distribution(0, 0xFFFFFFFFFFFFFFF);
	std::vector<long> int_vec(num_elements);
	for (int i = 0; i < num_elements; i++) {
		int_vec[i] = i;
	}
	std::random_shuffle(int_vec.begin(), int_vec.end());
	return int_vec;
}

/*
int main(int arg_num, char* args) {
	int positive_elements = 1000;
	int negative_elements = 1000000;
	double simple_fpr = 0;
	double amp_fpr = 0;
	double simple_LF = 0;
	double amp_LF = 0;
	int num_reps = 10;
	int bits_per = 8;
	int num_amplification = 1000000;
	for (int reps = 0; reps < num_reps; reps++) {
		std::vector<long> ints = generate_ints(positive_elements+negative_elements);
		BloomFilter best_filter = BloomFilter(0,0,0);
		double best_load_factor = 1;
		for (int i = 0; i < num_amplification; i++) {
			BloomFilter temp_filter = BloomFilter(positive_elements * bits_per, (double)bits_per * log(2), rand());
			for (int j = 0; j < positive_elements; j++) temp_filter.addElement(ints[j]);
			if (temp_filter.GetLoadFactor() < best_load_factor) {
				best_filter = temp_filter;
				best_load_factor = temp_filter.GetLoadFactor();
			}
		}
		BloomFilter simple_filter = BloomFilter(positive_elements * bits_per, (double)bits_per * log(2), rand());
		for (int i = 0; i < positive_elements; i++) simple_filter.addElement(ints[i]);
		size_t simple_fps = 0;
		size_t amp_fps = 0;
		for (int i = positive_elements; i < positive_elements + negative_elements; i++) if (simple_filter.testElement(ints[i]) == 1) simple_fps++;
		for (int i = positive_elements; i < positive_elements + negative_elements; i++) if (best_filter.testElement(ints[i]) == 1) amp_fps++;
		simple_fpr += (double)(simple_fps) / (double)(negative_elements);
		amp_fpr += (double)(amp_fps) / (double)(negative_elements);
		simple_LF += simple_filter.GetLoadFactor();
		amp_LF += best_filter.GetLoadFactor();
		printf("Trial Amp FPR:%f Trial Simple FPR:%f  Amp LF:%f Simple LF:%f\n", (double)(amp_fps) / (double)(negative_elements), (double)(simple_fps) / (double)(negative_elements), best_filter.GetLoadFactor(), simple_filter.GetLoadFactor());
	}
	simple_fpr = simple_fpr / num_reps;
	amp_fpr = amp_fpr / num_reps;
	simple_LF = simple_LF / num_reps;
	amp_LF = amp_LF / num_reps;
	printf("Num Positive:%ld Num Negative:%ld simple_fpr:%f amp_fpr:%f simple_LF:%f amp_LF:%f\n", positive_elements, negative_elements, simple_fpr, amp_fpr, simple_LF, amp_LF);
}*/


size_t getHash1(long x, int seed) {
	return CityHash64WithSeed((char*)&x, 4, 1234567+seed);
}

size_t getHash2(long x, int seed) {
	return CityHash64WithSeed((char*)&x, 4, 7654321+seed);
}

// Double hashing strategy recommended in Mitzenmacher Paper.
size_t getNthHash(size_t hash1, size_t hash2, size_t hash_num, size_t filter_size) {
	return (hash1 + hash_num * hash2) % filter_size;
}

int main(int arg_num, char* args) {
	int num_positives = 1000000;
	int num_negatives = 10000000;
	double simple_fpr = 0;
	double amp_fpr = 0;
	double simple_LF = 0;
	double amp_LF = 0;
	double prob_full_solo = 0;
	int num_reps = 10;
	double bits_per = 6;
	double bits_per_backup = 14;
	int num_hashes = (int)((double)bits_per*log(2))+3;
	int non_collisions = (int)((double)bits_per*log(2))-1;
	for (int outer = 0; outer < 20; outer++) {
		double standard_fpr_total = 0;
		double removed_fpr_total = 0;
		for (int reps = 0; reps < num_reps; reps++) {
			int seed = rand();
			std::vector<unsigned int> filter_(num_positives*bits_per);
			std::vector<long> ints = generate_ints(num_positives + num_negatives);
			for (int i = 0; i < num_positives; i++) {
				size_t hash1 = getHash1(ints[i], seed);
				size_t hash2 = getHash2(ints[i], seed);
				for (int j = 0; j < num_hashes; j++)filter_[getNthHash(hash1, hash2, j, filter_.size())]++;
			}
			std::vector<unsigned int> num_solo(num_positives, 0);
			for (int i = 0; i < num_positives; i++) {
				size_t hash1 = getHash1(ints[i], seed);
				size_t hash2 = getHash2(ints[i], seed);
				for (int j = 0; j < num_hashes; j++)num_solo[i] += (filter_[getNthHash(hash1, hash2, j, filter_.size())] == 1);
			}
			int full_solo = 0;
			for (int i = 0; i < num_positives; i++)if (num_solo[i] >= non_collisions)full_solo++;
			printf("Full Solo Prob:%f ", (double)full_solo / num_positives);
			std::vector<bool> filter_removed(num_positives*bits_per, 0);
			std::vector<bool> filter_removed_backup(full_solo*bits_per_backup, 0);
			std::vector<bool> filter_standard(num_positives*bits_per + bits_per_backup * full_solo, 0);
			int num_hashes_standard = (double)(num_positives*bits_per + bits_per_backup * full_solo) / num_positives * log(2);
			int num_hashes_backup = (double)(bits_per_backup)* log(2);
			for (int i = 0; i < num_positives; i++) {
				size_t hash1 = getHash1(ints[i], seed);
				size_t hash2 = getHash2(ints[i], seed);
				for (int j = 0; j < num_hashes_standard; j++)filter_standard[getNthHash(hash1, hash2, j, filter_standard.size())] = 1;
				if (num_solo[i] < non_collisions) {
					for (int j = 0; j < num_hashes; j++) filter_removed[getNthHash(hash1, hash2, j, filter_removed.size())] = 1;
				}
				else {
					for (int j = 0; j < num_hashes_backup; j++) filter_removed_backup[getNthHash(hash1, hash2, j, filter_removed_backup.size())] = 1;
				}
			}
			double standard_fpr = 0;
			double removed_fpr = 0;
			double removed_front_fpr = 0;
			double removed_backup_fpr = 0;
			for (int i = 0; i < num_negatives; i++) {
				size_t hash1 = getHash1(ints[num_positives + i], seed);
				size_t hash2 = getHash2(ints[num_positives + i], seed);
				bool removed_fps = true;
				bool backup_fps = true;
				bool standard_fps = true;
				for (int j = 0; j < num_hashes_standard; j++) {
					if (filter_standard[getNthHash(hash1, hash2, j, filter_standard.size())] == 0) {
						standard_fps = false;
						break;
					}
				}
				for (int j = 0; j < num_hashes; j++) {
					if (filter_removed[getNthHash(hash1, hash2, j, filter_removed.size())] == 0) {
						removed_fps = false;
						break;
					}
				}
				for (int j = 0; j < num_hashes_backup; j++) {
					if (filter_removed_backup[getNthHash(hash1, hash2, j, filter_removed_backup.size())] == 0) {
						backup_fps = false;
						break;
					}
				}
				if (backup_fps == true)removed_backup_fpr++;
				if (removed_fps == true)removed_front_fpr++;
				if (backup_fps == true || removed_fps == true) removed_fpr++;
				if (standard_fps == true)standard_fpr++;
			}
			standard_fpr /= num_negatives;
			removed_fpr /= num_negatives;
			removed_front_fpr /= num_negatives;
			removed_backup_fpr /= num_negatives;
			double standard_LF = 0;
			double removed_LF = 0;
			double backup_LF = 0;
			for (int i = 0; i < filter_standard.size(); i++)standard_LF += filter_standard[i];
			for (int i = 0; i < filter_removed.size(); i++)removed_LF += filter_removed[i];
			for (int i = 0; i < filter_removed_backup.size(); i++)backup_LF += filter_removed_backup[i];
			standard_LF /= filter_standard.size();
			removed_LF /= filter_removed.size();
			backup_LF /= filter_removed_backup.size();
			printf("TRIAL %d standard_fpr:%f removed_fpr:%f standard_LF:%f removed_LF:%f backup_LF:%f\n", reps, standard_fpr, removed_fpr, standard_LF, removed_LF, backup_LF);
			printf("removed_front_fpr:%f removed_backup_fpr:%f removed_front_size:%ld removed_backup_size:%ld\n", removed_front_fpr, removed_backup_fpr, filter_removed.size(), filter_removed_backup.size());
			standard_fpr_total += standard_fpr;
			removed_fpr_total += removed_fpr;
		}
		standard_fpr_total /= num_reps;
		removed_fpr_total /= num_reps;
		printf("standard_fpr:%f removed_fpr:%f\n", standard_fpr_total, removed_fpr_total);
	}
}