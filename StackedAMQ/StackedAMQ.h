#pragma once
#include <vector>

class StackedAMQ
{
private:
	std::vector<std::vector<std::vector<double>>> opt_table;
public:
	StackedAMQ(std::vector<std::pair<char*, unsigned int>> pos_elements, std::vector<std::pair<char*, unsigned int>> neg_elements, const char* AMQ_Type, size_t total_size, bool opt_table_exists, const char* opt_table_file_path);
	~StackedAMQ();
	bool LookupElement(char* element, unsigned int length);
	void InsertElement(char* element, unsigned int length);
	void DeleteElement(char* element, unsigned int length);
	void GenerateOptTable(const char* AMQ_Type, double granularity);
	void SaveOptTable(const char* filepath);
	void LoadOptTable(const char* filepath);
	void PrintOptTable();
	size_t GetSize();
	double GetFPR();
	std::vector<double> CalculateLayerFPRs();
	std::vector<double> RetrieveLayerFPRsFromOptTable();
};

