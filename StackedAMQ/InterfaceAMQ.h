#pragma once
#include <cstddef>

class InterfaceAMQ
{
public:
	double fpr_;
	size_t total_size_;
	size_t elements_added_;
	InterfaceAMQ(double fpr, size_t num_expected_elements);
	InterfaceAMQ(size_t total_size_);
	~InterfaceAMQ();
	virtual void InsertElement();
	virtual bool LookupElement();


	// These functions must be overridden!
	static double SizeFunction(double fpr, size_t num_expected_elements);
	static double FPRFunction(size_t total_size, size_t num_expected_elements);
};

