#include "../StackedBloomFilter/StackedBloomFilter.h"
#include <fstream>
#include <stdlib.h>
#include <random>

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


int main(int arg_num, char* args) {
	std::ofstream fs;
	fs.open("C:\\Users\\Owner\\Documents\\School\\CS 91R\\Stacked_Filter_Data.csv");
	fs << "Beta,Psi,Bits Available,Equal Fprs,Num Layers,Used Bits,Total FPR,Known FPR,UnknownFPR,Filter Checks For Positive,Filter Checks For Negative\n";
	int total_elements = 1000000;
	for (double beta = .9; beta < .91; beta += .04) {
		for (double psi = .8; psi < .81; psi += .04) {
			for (double bits = 7; bits < 14; bits += .5) {
				int num_positives = total_elements * beta;
				int num_known_negatives = total_elements * (1 - beta);
				int num_unknown_negatives = num_known_negatives;
				printf("num_known_negatives:%d\n", num_known_negatives);
				int total_size = num_positives * bits;
				for (int i = 0; i <= 2; i++) {
					bool equal_fprs = 0;
					int num_layers = 0;
					if (i == 1)num_layers = 1;
					if (i == 2) num_layers = 2;
					if (i == 3) num_layers = 3;

					printf("beta =%f, psi=%f, bits = %f, equal_fprs=%d, num_layers=%d\n", beta, psi, bits, equal_fprs, num_layers);
					double known_fpr = 0;
					double unknown_fpr = 0;
					double total_fpr = 0;
					double used_bits = 0;
					double checks_per_pos = 0;
					double checks_per_neg = 0;
					int num_reps = 3;
					for (int reps = 0; reps < num_reps; reps++) {
						std::vector<long> ints = generate_ints(num_positives + num_known_negatives + num_unknown_negatives);
						std::vector<long> positives = std::vector<long>(ints.begin(), ints.begin() + num_positives);
						std::vector<long> known_negatives = std::vector<long>(ints.begin() + num_positives, ints.begin() + num_positives + num_known_negatives);
						std::vector<long> unknown_negatives = std::vector<long>(ints.begin() + num_positives + num_known_negatives, ints.end());
						StackedBloomFilter filter = StackedBloomFilter(num_layers, positives, known_negatives, total_size, psi, .0000001, equal_fprs);
						filter.PrintLayerDiagnostics();
						filter.ResetNumFilterChecks();
						for (int i = 0; i < num_positives; i++) filter.TestElement(positives[i]);
						checks_per_pos += (double)filter.NumFilterChecks() / num_positives;
						filter.ResetNumFilterChecks();
						int kfp = 0;
						for (int i = 0; i < num_known_negatives; i++) {
							if (filter.TestElement(known_negatives[i])) {
								kfp++;
							}
						}
						checks_per_neg += (double)filter.NumFilterChecks() / num_known_negatives * psi;
						filter.ResetNumFilterChecks();
						int ukfp = 0;
						for (int i = 0; i < num_unknown_negatives; i++) {
							if (filter.TestElement(unknown_negatives[i])) {
								ukfp++;
							}
						}
						checks_per_neg += (double)filter.NumFilterChecks() / num_unknown_negatives * (1 - psi);
						known_fpr += (double)(kfp) / (double)(num_known_negatives);
						unknown_fpr += (double)(ukfp) / (double)(num_unknown_negatives);
						total_fpr += psi * (double)(kfp) / (double)(num_known_negatives)+ (1 - psi)*(double)(ukfp) / (double)(num_unknown_negatives);
						used_bits += (double)filter.TotalSize() / num_positives;
						printf("Trial FPR:%f\n", psi * (double)(kfp) / (double)(num_known_negatives)+(1 - psi)*(double)(ukfp) / (double)(num_unknown_negatives));
					}
					known_fpr = known_fpr / num_reps;
					unknown_fpr = unknown_fpr / num_reps;
					total_fpr = total_fpr / num_reps;
					used_bits = used_bits / num_reps;
					checks_per_neg = checks_per_neg / num_reps;
					checks_per_pos = checks_per_pos / num_reps;
					fs << beta << "," << psi << "," << bits << "," << equal_fprs << "," << num_layers << ","<< used_bits << "," << total_fpr << "," << known_fpr << "," << unknown_fpr << ","<<checks_per_pos<<","<<checks_per_neg<<"\n";
					printf(", fpr=%f, pos_checks=%f, neg_checks=%f", total_fpr, checks_per_pos, checks_per_neg);
					printf(" Used Bits= %f\n", used_bits);
				}
			}
		}
	}
	fs.close();
}