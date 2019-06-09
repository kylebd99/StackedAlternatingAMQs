#pragma once
#include "Common.h"
#include "BloomFilter.h"
#include "C:\Users\Owner\Documents\School\CS 91R\nlopt\nlopt.h"

class StackedBloomFilter
{
public:
	unsigned int num_layers_;
	double fpr_;
	double beta_;
	double psi_;
	double penalty_coef_;
	size_t num_positive_;
	size_t num_negative_;
	size_t total_size_;
	bool equal_layer_fprs_ = false;
	std::vector<BloomFilter> filter_array_;
	std::vector<double> filter_fprs_;
	StackedBloomFilter(int num_filters, std::vector<long> positives, std::vector<long> negatives, size_t total_size, double psi, double penalty_coef, bool equal_layer_fprs, std::vector<double> layer_fprs);
	StackedBloomFilter(int num_filters, std::vector<long> positives, std::vector<long> negatives, size_t total_size, double psi, double penalty_coef, bool equal_layer_fprs) : StackedBloomFilter(num_filters, positives, negatives, total_size, psi, penalty_coef, equal_layer_fprs, std::vector<double>()) {}
	StackedBloomFilter(int num_filters, std::vector<long> positives, std::vector<long> negatives, size_t total_size, double psi, double penalty_coef) : StackedBloomFilter(num_filters, positives, negatives, total_size, psi, penalty_coef, false) {}
	StackedBloomFilter();
	std::vector<double> EstimateLayerFPR();
	void AddPositiveElement(int element);
	bool TestElement(int element);
	size_t TotalSize();
	size_t NumFilterChecks();
	void ResetNumFilterChecks();
	void PrintLayerDiagnostics();

	static double FprFunctionVaried(unsigned num_filters, const double* lfprs, double* grad, void* filterptr) {
		StackedBloomFilter* filter = (StackedBloomFilter*)filterptr;
		double beta = filter->beta_;
		size_t total_size = filter->total_size_;
		int num_layers = filter->num_layers_;
		double psi = filter->psi_;
		double known_fpr = lfprs[0];
		double penalty_coef = filter->penalty_coef_;
		for (int i = 1; i <= num_layers; i++) {
			known_fpr = known_fpr * lfprs[2 * i];
		}
		double unknown_fpr_side = 0;
		for (int i = 1; i <= num_layers; i++) {
			double temp_fpr = lfprs[0];
			for (int j = 1; j <= 2 * (i - 1); j++) {
				temp_fpr = temp_fpr * lfprs[j];
			}
			unknown_fpr_side += temp_fpr * (1 - lfprs[2 * i - 1]);
		}
		double unknown_fpr_end = lfprs[0];
		for (int i = 1; i <= 2 * num_layers; i++) {
			unknown_fpr_end = unknown_fpr_end * lfprs[i];
		}
		double total_fpr = psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end);
		// Penalty Function For Number of Hashes
		int num_hashes = 0;
		for (int i = 0; i<num_layers * 2; i++) num_hashes += -log(lfprs[i]) / log(2);
		if (grad != NULL) {
			grad[0] = (psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end))/lfprs[0];
			for (int k = 1; k < 2 * num_layers + 1; k++) {
				grad[k] = 0;
				if (k % 2 == 0) {
					double temp_fpr = 1;
					for (int i = 0; i <= num_layers; i++) temp_fpr = temp_fpr * lfprs[2 * i];
					grad[k] += psi * temp_fpr / lfprs[k];
					for (int i = 1; i <= num_layers; i++) {
						if (2 * (i - 1) >= k) {
							temp_fpr = 1;
							for (int j = 0; j <= 2 * (i - 1); j++) {
								temp_fpr = temp_fpr * lfprs[j];
							}
							grad[k] += (1 - psi)*temp_fpr*(1 - lfprs[2 * i - 1]) / lfprs[k];
						}
					}
					temp_fpr = 1;
					for (int i = 0; i < num_layers * 2 + 1; i++) temp_fpr = temp_fpr * lfprs[i];
					grad[k] += temp_fpr / lfprs[k];
				}
				else {
					double temp_fpr = 1;
					for (int i = 1; i <= num_layers; i++) {
						if (2 * i - 1 == k) {
							temp_fpr = 1;
							for (int j = 0; j <= 2 * (i - 1); j++) temp_fpr = temp_fpr * lfprs[j];
							grad[k] += -temp_fpr;
						}
						else if (2 * i - 1 > k) {
							temp_fpr = 1;
							for (int j = 0; j <= 2 * (i - 1); j++) temp_fpr = temp_fpr * lfprs[j];
							grad[k] += temp_fpr * (1 - lfprs[2 * i - 1]) / lfprs[k];
						}
					}
					temp_fpr = 1;
					for (int j = 0; j < 2 * num_layers + 1; j++) temp_fpr = temp_fpr * lfprs[j];
					grad[k] += temp_fpr / lfprs[k];
					grad[k] = grad[k] * (1 - psi);
				}
			}
			// Adding in the penalty function to the Gradient
			for (int k = 0; k < num_layers * 2 + 1; k++) grad[k] = grad[k] * (1 + num_hashes * penalty_coef + total_fpr * (-1 / lfprs[k] / log(2)*penalty_coef));
		}
//		printf("lfprs[0]:%f lfprs[1]:%f lfprs[2]:%f fpr:%f\n", lfprs[0], lfprs[1], lfprs[2], psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end));
		return (psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end))*(1+num_hashes*penalty_coef);
	};

	static double FprFunctionEqual(unsigned num_filters, const double* lfprs, double* grad, void* filterptr) {
		StackedBloomFilter* filter = (StackedBloomFilter*)filterptr;
		double beta = filter->beta_;
		int num_layers = filter->num_layers_;
		double psi = filter->psi_;
		double known_fpr = lfprs[0];
		double penalty_coef = filter->penalty_coef_;
		for (int i = 1; i <= num_layers; i++) {
			known_fpr = known_fpr * lfprs[0];
		}
		double unknown_fpr_side = 0;
		for (int i = 1; i <= num_layers; i++) {
			double temp_fpr = lfprs[0];
			for (int j = 1; j <= 2 * (i - 1); j++) {
				temp_fpr = temp_fpr * lfprs[0];
			}
			unknown_fpr_side += temp_fpr * (1 - lfprs[0]);
		}
		double unknown_fpr_end = lfprs[0];
		for (int i = 1; i <= 2 * num_layers; i++) {
			unknown_fpr_end = unknown_fpr_end * lfprs[0];
		}
		int num_hashes = 0;
		for (int i = 0; i<num_layers * 2; i++) num_hashes += -log(lfprs[i]) / log(2);
		/*if (grad != NULL) {
			grad[0] += psi * (num_layers + 1)*pow(lfprs[0], num_layers);
			double temp_fpr = 0;
			for (int i = 1; i <= num_layers; i++) grad[0] += (1-psi)*((2 * i - 1)*pow(lfprs[0], 2 * (i - 1)) - 2 * i*pow(lfprs[0], 2 * i - 1));
			grad[0] += (1 - psi)*((2 * num_layers + 1)*pow(lfprs[0], 2 * num_layers));
		}*/
		//printf("lfprs[0]:%f lfprs[1]:%f lfprs[2]:%f fpr:%f\n", lfprs[0], lfprs[0], lfprs[0], psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end));
		return (psi * known_fpr + (1 - psi)*(unknown_fpr_side + unknown_fpr_end))*(1 + num_hashes * penalty_coef);
	};

	static double SizeFunctionVaried(unsigned num_filters, const double* lfprs, double* grad, void* filterptr) {
		StackedBloomFilter* filter = (StackedBloomFilter*)filterptr;
		double beta = filter->beta_;
		size_t total_size = filter->total_size_*.998;
		int num_layers = filter->num_layers_;
		size_t num_positive = filter->num_positive_;
		size_t num_negative = filter->num_negative_;
		double size = 0;
		/*for (int i = 1; i <= num_layers; i++){
			double negative_fpr = lfprs[0];
			for (int j = 1; j <= i - 1; j++) {
				negative_fpr = negative_fpr * lfprs[2 * j];
			}
			size += negative_fpr * log(lfprs[2 * i - 1]) * (1 - beta);
			double positive_fpr = lfprs[1];
			for (int j = 1; j <= i - 1; j++) {
				positive_fpr = positive_fpr * lfprs[2 * j + 1];
			}
			size += positive_fpr * log(lfprs[2 * i]) * beta;
		}
		size = -size * (num_negative + num_positive) / log(2) / log(2);*/
		double positive_fpr = 1;
		double negative_fpr = 1;
		for (int i = 0; i <= num_layers*2; i++) {
			int num_hashes = max((int)(round(-log(lfprs[i]) / log(2)) + .5),1);
			double temp_size;
			if (((i+2) % 2) == 0) {
				temp_size = 1 / (1 - (double)pow((1 - (double)pow(lfprs[i], (double)1 / num_hashes)), (double)1 / (positive_fpr*num_positive*num_hashes)));
				negative_fpr *= lfprs[i];
			}
			else {
				temp_size = 1 / (1 - (double)pow((1 - (double)pow(lfprs[i],(double) 1 / num_hashes)), (double)1 / (negative_fpr*num_negative*num_hashes)));
				positive_fpr *= lfprs[i];
			}
			// Bloom Filter FPR formulas do not work well on small filters, so we put a floor on the size of each layer.
			if (temp_size < 2000 && lfprs[i]<.99) temp_size = 2000;
			size += temp_size;
		}/*
		if (grad != NULL) {
			// Calculating the gradients for the various layer fprs.
			size_t num_elements = num_positive + num_negative;
			double l2 = 1 / log(2) / log(2);
			for (int k =  0; k < num_layers * 2 + 1; k++) {
				grad[k] = 0;
				if (k % 2 == 0) {
					double temp_fpr = 1;
					for (int j = 1; j <= k / 2; j++) {
						temp_fpr = temp_fpr * lfprs[2 * j - 1];
					}
					grad[k] += -beta * temp_fpr*num_elements / lfprs[k] * l2;
					for (int i = k / 2 + 1; i <= num_layers; i++) {
						temp_fpr = 1;
						for (int j = 0; j <= i - 1; j++) {
							temp_fpr = temp_fpr * lfprs[2 * j];
						}
						grad[k] += -(1 - beta)*temp_fpr*log(lfprs[2 * i - 1])*num_elements / lfprs[k] * l2;
					}
				}
				else {
					double temp_fpr = 1;
					for (int j = 1; j <= (k + 1) / 2 - 1; j++) {
						temp_fpr = temp_fpr * lfprs[2 * j];
					}
					grad[k] += -(1 - beta)*temp_fpr*num_elements / lfprs[k] * l2;
					for (int i = (k + 1) / 2; i <= num_layers; i++) {
						temp_fpr = 1;
						for (int j = 1; j <= i; j++) {
							temp_fpr = temp_fpr * lfprs[2 * j - 1];
						}
						grad[k] += -beta * temp_fpr*log(lfprs[i])*num_elements / lfprs[k] * l2;
					}
				}
			}
		}*/
//		printf("Size Function:%d FPR0:%f FPR1:%f FPR2:%f\n", size - total_size, lfprs[0], lfprs[1], lfprs[2]);
		return size - total_size;
	};

	static double SizeFunctionEqual(unsigned num_filters, const double* lfprs, double* grad, void* filterptr) {
		StackedBloomFilter* filter = (StackedBloomFilter*)filterptr;
		double beta = filter->beta_;
		size_t total_size = filter->total_size_*.999;
		int num_layers = filter->num_layers_;
		size_t num_positive = filter->num_positive_;
		size_t num_negative = filter->num_negative_;
		double size = 0;
		double positive_fpr = 1;
		double negative_fpr = 1;
		for (int i = 0; i <= num_layers*2; i++) {
			int num_hashes = max((int)(round(-log(lfprs[i]) / log(2))), 1);
			if ((i % 2) == 0) {
				size += 1 / (1 - (double)pow((1 - (double)pow(lfprs[i], (double)1 / num_hashes)), (double)1 / (positive_fpr*num_positive*num_hashes)));
				negative_fpr *= lfprs[i];
			}
			else {
				size += 1 / (1 - (double)pow((1 - (double)pow(lfprs[i], (double)1 / num_hashes)), (double)1 / (negative_fpr*num_negative*num_hashes)));
				positive_fpr *= lfprs[i];
			}
		}/*if (grad != NULL) {
			// Calculating the gradients for the various layer fprs.
			size_t num_elements = num_positive + num_negative;
			double l2 = 1 / log(2) / log(2);
			double temp_fpr = 0;
			for (int i = 1; i <= num_layers; i++) temp_fpr += pow(lfprs[0], i-1) * i;
			grad[0] += -(beta + temp_fpr)*log(lfprs[0])*num_elements*l2;
			temp_fpr = 0;
			for (int i = 1; i <= num_layers; i++) temp_fpr += pow(lfprs[0], i);
			grad[0] += -(beta + temp_fpr) / lfprs[0] * num_elements*l2;
		}*/
		//printf("Size Function:%d FPR:%f\n", size - total_size + 1, lfprs[0]);
		return size - total_size;
	};

};

