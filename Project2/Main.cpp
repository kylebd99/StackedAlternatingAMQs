#include "../StackedBloomFilter/StackedBloomFilter.h"



void main(int arg_num, char* args) {
	StackedBloomFilter filter;
	filter.penalty_coef_ = .001;
	filter.beta_ = .05;
	filter.psi_ = .65;
	filter.num_layers_ = 1;
	filter.num_negative_ = 35000;
	filter.num_positive_ = 65000;
	filter.total_size_ = filter.num_positive_ * 14;
	filter.filter_fprs_ = std::vector<double>(3, 0.00123488245160775800);
	std::vector<double> grads(3,0);
	double size = filter.SizeFunctionVaried(3, filter.filter_fprs_.data(), grads.data(), &filter);
	double fpr = filter.FprFunctionVaried(3, filter.filter_fprs_.data(), grads.data(), &filter);
	printf("FPR = %f , Size = %f, grad0 = %f, grad1 = %f, grad2 = %f\n", fpr, size, grads[0], grads[1], grads[2]);
	int x;
}

