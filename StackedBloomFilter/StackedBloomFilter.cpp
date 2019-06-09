#include "StackedBloomFilter.h"

StackedBloomFilter::StackedBloomFilter() {
}

StackedBloomFilter::StackedBloomFilter(int num_layers, std::vector<long> positives, std::vector<long> negatives, size_t total_size, double psi, double penalty_coef, bool equal_layer_fprs, std::vector<double> layer_fpr)
{
	penalty_coef_ = penalty_coef;
	equal_layer_fprs_ = equal_layer_fprs;
	num_layers_ = num_layers;
	num_positive_ = positives.size();
	num_negative_ = negatives.size();
	total_size_ = total_size;
	psi_ = psi;
	beta_ = (double)num_positive_ / (double)(num_positive_ + num_negative_);
	if(layer_fpr.size()==0) layer_fpr = EstimateLayerFPR();
	filter_array_ = std::vector<BloomFilter>();
	// Preallocate each of the bloom filters to their estimated sizes
	double positive_fpr = 1;
	double negative_fpr = 1;
	for (int i = 0; i <= num_layers * 2; i++) {
		int num_hashes = max((int)(round(-log(layer_fpr[i]) / log(2))), 1);
		size_t size = 0;
		if ((i % 2) == 0) {
			size = 1 / (1 - (double)pow((1 - (double)pow(layer_fpr[i], (double)1 / num_hashes)), (double)1 / (positive_fpr*num_positive_*num_hashes)));
			negative_fpr *= layer_fpr[i];
			}
		else {
			size = 1 / (1 - (double)pow((1 - (double)pow(layer_fpr[i], (double)1 / num_hashes)), (double)1 / (negative_fpr*num_negative_*num_hashes)));
			positive_fpr *= layer_fpr[i];
		}
		// Bloom Filter FPR formulas do not work well on small filters, so we put a floor on the size of each layer.
		if (size < 2000 && layer_fpr[i]<.99) size = 2000;
		filter_array_.push_back(BloomFilter(size, num_hashes, rand()));
	}
	for (int i = 0; i < num_positive_; i++) {
		filter_array_[0].addElement(positives[i]);
	}
	// Add the elements to the filters in a descending order by layer (filter-pair)
	size_t num_negative_fp = num_negative_;
	size_t num_positive_fp = num_positive_;
	std::vector<int> negative_fp(negatives.begin(), negatives.end());
	std::vector<int> positive_fp(positives.begin(), positives.end());
	for (int i = 0; i < num_layers_; i++) {
		int temp_neg_fp = 0;
		for (int j = 0; j < num_negative_fp; j++) {
			if (filter_array_[2 * i].testElement(negative_fp[j])) {
				negative_fp[temp_neg_fp] = negative_fp[j];
				temp_neg_fp++;
			}
		}
		num_negative_fp = temp_neg_fp;
		for (int j = 0; j < num_negative_fp; j++) {
			filter_array_[2 * i+1].addElement(negative_fp[j]);
		}
		int temp_pos_fp = 0;
		for (int j = 0; j < num_positive_fp; j++) {
			if (filter_array_[2 * i+1].testElement(positive_fp[j])) {
				positive_fp[temp_pos_fp] = positive_fp[j];
				temp_pos_fp++;
			}
		}
		num_positive_fp = temp_pos_fp;
		for (int j = 0; j < num_positive_fp; j++) {
			filter_array_[2 * i + 2].addElement(positive_fp[j]);
		}
	}
}

std::vector<double> StackedBloomFilter::EstimateLayerFPR() {
	double one_level_fpr = exp(-(long long)total_size_*log(2)*log(2) / beta_ / (num_negative_ + num_positive_));
	int num_fprs = num_layers_ * 2 + 1;
	std::vector<double> one_level_fprs = std::vector<double>(num_fprs, one_level_fpr);
	for (int i = 1; i < num_fprs; i++) one_level_fprs[i] = 1;
	if (num_layers_ == 0) {
		filter_fprs_ = one_level_fprs;
		return filter_fprs_;
	}
	double* zeros = (double*)calloc(num_fprs, sizeof(double));
	for (int i = 0; i < num_fprs; i++) zeros[i] = 0.00000000000000000001;
	double* ones = (double*)calloc((num_fprs), sizeof(double));
	for (int i = 0; i < num_fprs; i++) ones[i] = 1;
	nlopt_opt equal_fpr_opt = nlopt_create(NLOPT_GN_ISRES, 1);
	nlopt_set_lower_bounds(equal_fpr_opt, zeros);
	nlopt_set_upper_bounds(equal_fpr_opt, ones);
	nlopt_set_maxtime(equal_fpr_opt, .25);
	nlopt_add_inequality_constraint(equal_fpr_opt, &StackedBloomFilter::SizeFunctionEqual, this, total_size_*.0005);
	nlopt_set_min_objective(equal_fpr_opt, &StackedBloomFilter::FprFunctionEqual, this);
	double equal_fpr_fpr = 0;
	std::vector<double> lfprs(1, .5);
	if (nlopt_optimize(equal_fpr_opt, lfprs.data(), &equal_fpr_fpr) < 0) printf("Equal Opt Error!!\n");
	std::vector<double> equal_fprs = std::vector<double>(num_layers_ * 2 + 1, lfprs[0]); 
	if (equal_layer_fprs_) {
		filter_fprs_ = equal_fprs;
		return filter_fprs_;
	}
	if (equal_fpr_fpr > one_level_fpr) {
		filter_fprs_ = one_level_fprs;
		printf("Using One Level Start\n");
	}
	else {
		printf("Using Equal Level Start\n");
		filter_fprs_ = equal_fprs;
	}
	nlopt_opt global_fpr_opt = nlopt_create(NLOPT_GN_ISRES, num_fprs);
	nlopt_set_lower_bounds(global_fpr_opt, zeros);
	nlopt_set_upper_bounds(global_fpr_opt, ones);
	nlopt_set_maxtime(global_fpr_opt, .5 + 1*(num_layers_ - 1)); //Providing more time for more parameters
	nlopt_add_inequality_constraint(global_fpr_opt, &StackedBloomFilter::SizeFunctionVaried, this, total_size_*.0005);
	nlopt_set_min_objective(global_fpr_opt, &StackedBloomFilter::FprFunctionVaried, this);
	double variable_fpr_fpr = 1;
	nlopt_result global_ret_status = nlopt_optimize(global_fpr_opt, filter_fprs_.data(), &variable_fpr_fpr);
	if (global_ret_status == -4) printf("ERROR!!!!!  Roundoff Errors Reached in Global Optimization\n");
	else if (global_ret_status < 0) printf("ERROR!!! GENERALLY SPEAKING IN GLOBAL\n");
	nlopt_opt local_fpr_opt = nlopt_create(NLOPT_LN_COBYLA, num_fprs);
	nlopt_set_lower_bounds(local_fpr_opt, zeros);
	nlopt_set_upper_bounds(local_fpr_opt, ones);
	nlopt_set_maxtime(local_fpr_opt, .5 + 2*(num_layers_ - 1)); //Providing more time for more parameters
	nlopt_set_ftol_rel(local_fpr_opt, .00001);
	nlopt_add_inequality_constraint(local_fpr_opt, &StackedBloomFilter::SizeFunctionVaried, this, total_size_*.0005);
	nlopt_set_min_objective(local_fpr_opt, &StackedBloomFilter::FprFunctionVaried, this);
	variable_fpr_fpr = 1;
	nlopt_result local_ret_status = nlopt_optimize(local_fpr_opt, filter_fprs_.data(), &variable_fpr_fpr);
	if (local_ret_status == -4) printf("ERROR!!!!!  Roundoff Errors Reached in Local Optimization\n");
	else if (local_ret_status < 0) printf("ERROR!!! GENERALLY SPEAKING IN LOCAL\n");
	return filter_fprs_;
}

bool StackedBloomFilter::TestElement(int element) {
	if (filter_array_[0].testElement(element) == false) {
		return false;
	}
	for (int i = 0; i < num_layers_; i++) {
		if (filter_array_[2 * i + 1].testElement(element) == false) return true;
		if (filter_array_[2 * i + 2].testElement(element) == false) return false;
	}
	return true;
}

void StackedBloomFilter::AddPositiveElement(int element) {
	filter_array_[0].addElement(element);
	for (int i = 0; i < num_layers_; i++) {
		if (filter_array_[2 * i + 1].testElement(element) == true) {
			filter_array_[2 * i + 2].addElement(element);
		}
		else {
			return;
		}
	}
}

size_t StackedBloomFilter::TotalSize() {
	size_t size = 0;
	size += filter_array_[0].GetSize();
	for (int i = 0; i < num_layers_; i++) {
		size += filter_array_[2 * i + 1].GetSize();
		size += filter_array_[2 * i + 2].GetSize();
	}
	return size;
}

size_t StackedBloomFilter::NumFilterChecks() {
	size_t num_filter_checks = 0;
	for (int i = 0; i < num_layers_ * 2 + 1; i++) {
		num_filter_checks += filter_array_[i].num_checks_;
	}
	return num_filter_checks;
}

void StackedBloomFilter::ResetNumFilterChecks() {
	for (int i = 0; i < num_layers_ * 2 + 1; i++) filter_array_[i].num_checks_ = 0;
}

void StackedBloomFilter::PrintLayerDiagnostics() {
	for (int i = 0; i < num_layers_ * 2 + 1; i++) {
		printf("Layer %d FPR: %.20f Size:%d Num Hashes:%d Num Elements:%ld Load Factor:%f \n", i, filter_fprs_[i], filter_array_[i].GetSize(), filter_array_[i].num_hashes_, filter_array_[i].GetNumElements(), filter_array_[i].GetLoadFactor());
	}
}
