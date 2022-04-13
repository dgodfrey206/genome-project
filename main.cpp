#include "UnorderedMap.hpp"
#include "blast.hpp"
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>
#include <iterator>
#include <random>
#include <unordered_map>
#include <chrono>
#include <iostream>
#include <cassert>
#include <sstream>
#include <utility>
#include <cstring>
using namespace std;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;

static const int WORD_SIZE = Blast_DB::WORD_SIZE;
static const int SENTENCE_SIZE = 50;

template<typename T>
T roundFloorMultiple( T value, T multiple )
{
    if (multiple == 0) return value;
    return static_cast<T>(std::floor(static_cast<double>(value)/static_cast<double>(multiple))*static_cast<double>(multiple));
}

using Data = Blast_DB::data;

/*
Go through each 50-mer string in sample_hw_dataset.txt, break those strings down into 11-mer 
strings, compare them with the hash table and store the matches in a stack. For each match that 
we found, use seed expansion to get the full 50-mer string from the test_genome.txt file. Using 
these two 50-mer strings, send them to the Needleman Wunsch algorithm, and output the result.
*/
using namespace std;
void ProcessDataset(std::string genome, std::string file, int iterations) {
	Blast_DB db(genome);
	db.store_polymers();
	auto& table = db.table();
	std::cout << "1c " << iterations << "\n";
	std::ifstream test(file);
	assert(test.is_open());
	std::vector<Data> stk;
	int pHits = 0;
	UnorderedMapPool found;
	for (std::string str; std::getline(test, str); ) {
		if (str[0] != '>') {
			for (std::size_t i = 0; i < str.size(); i++) {
				std::string word = str.substr(i, WORD_SIZE);
				if (found[word] == 0 && word.size() == WORD_SIZE && table.find(word) != table.end()) {
					found[word] = 1;
					stk.push_back({ str, i, table[word] });
					int idx = table[word] - i;
					if (idx >= 0) {
						std::string genome_substr = genome.substr(idx, 50);
						std::cout << genome_substr << " " << str << '\n';
						auto p = Blast_DB::query(genome_substr, str);
						std::cout << "Genome location for best hit: " << idx << '\n';
						std::cout << "Score: " << p.first << '\n';
						if (p.first == 100) pHits++;
						auto p2 = p.second;
						std::cout << p2.first << '\n';
						for (int i = 0; i < p2.first.size(); i++) {
							if (p2.first[i] != '-' && p2.second[i] != '-') {
								if (p2.first[i] != p2.second[i])
									std::cout << "x";
								else
									std::cout << "|";
							} else {
								std::cout << " ";
							}
						}
						std::cout << '\n';
						std::cout << p2.second << "\n\n";
					}
				}
			}
		}
	}
	std::cout << "Perfect hits: " << pHits << '\n';
}

template<class T>
void print(std::vector<std::vector<T>> arr) {
  std::cout << "[\n";
  
  for (int i=0; i<arr.size(); i++) {
    const char* sep = "";
    for (int j=0; j<arr[i].size(); j++) {
      std::cout << sep << arr[i][j];
      sep = ", ";
    }
  }
  std::cout<<"]\n";
}

void runq1(int iterations, std::string genome, std::vector<Data>& stk) {
	Blast_DB db(genome.substr(0, iterations));
	std::cout << "Number of characters in the genome: " << genome.size() << '\n';
	assert(genome.size());
	auto& table = db.table();
	std::cout << "Number of " << WORD_SIZE << " character fragments possible: " << (genome.size() - WORD_SIZE + 1) << "\n";

	std::vector<int> q;
	std::cout << "Generating random queries...\n";

	std::cout << "Populating hash table...\n";
	db.store_polymers();
	std::cout << "Hash table populated\n";
	std::cout << "Starting timer:\n";

	std::cout << "Processing queries...\n";
	int count = 0;
	auto t1 = high_resolution_clock::now();
	
	UnorderedMapPool found;
	string last;
	std::string sentence = genome.substr(0, iterations);
	for (std::size_t i = 0; i < sentence.size(); i++) {
		std::string word = sentence.substr(i, WORD_SIZE);
		
		if (word.size() != WORD_SIZE) break;
		
		if (found[word] == 0 && table.find(word) != table.end()) {
			found[word] = 1;
			stk.push_back(Data{ word, i, table[word] });
			count++;
		}
	}
	
	auto t2 = high_resolution_clock::now();
	std::chrono::duration<double> ms_double = t2 - t1;
	std::cout << "Ended timer\n";
	std::cout << "Time taken: " << ms_double.count() << " sec\n";

	std::cout << "Total 11 fragments matched: " << count << "\n";
	std::cout << "Total queries used: " << iterations << "\n";
}

void runq2(int c, std::string genome, std::vector<Data>& stk) {
	Blast_DB db(genome.substr(0, c));
	std::cout << "Number of characters in the genome: " << genome.size() << '\n';
	assert(genome.size());
	auto& table = db.table();
	std::cout << "Number of " << WORD_SIZE << " character fragments possible: " << (genome.size() - WORD_SIZE + 1) << "\n";

	std::vector<int> q;
	std::cout << "Generating random queries...\n";

	std::random_device rd;
	std::mt19937 gen(rd());
	std::default_random_engine dre(rd());
	std::geometric_distribution<> d(0.05);

	for (int i = 0; i < c; i++)
		q.push_back(d(gen) % c);

	std::cout << "Populating hash table...\n";
	db.store_polymers();
	std::cout << "Hash table populated\n";
	std::cout << "Starting timer:\n";

	std::cout << "Processing queries...\n";
	int count = 0;
	int idx = 0;
	auto t1 = high_resolution_clock::now();
	UnorderedMapPool found;
	
	for (int i = 0; i < q.size(); i++) {
		idx += q[i];
		int newIdx = roundFloorMultiple(idx % c, 50);
		std::string sentence = genome.substr(newIdx, 50);
		if (sentence.size() != 50) continue;
		for (std::size_t i = 0; i < sentence.size(); i++) {
			std::string word = sentence.substr(i, WORD_SIZE);
			if (word.size() == WORD_SIZE) {
				if (found[word] == 0 && table.find(word) != table.end()) {
					found[word] = 1;
					stk.push_back(Data{ word, i, table[word] });
					count++;
				}
			}
		}
	}
	auto t2 = high_resolution_clock::now();
	std::chrono::duration<double> ms_double = t2 - t1;
	std::cout << "Ended timer\n";
	std::cout << "Time taken: " << ms_double.count() << " sec\n";

	std::cout << "Total 11 fragments matched: " << count << "\n";
	std::cout << "Total queries used: " << q.size() << "\n";
	std::cout << "Perfect hits(score = 100): " << '\n';
}

void q1(int c, std::string genome, std::vector<Data>& stk) {
	for (int i = 1; i <= 3; i++) {
		string x(i, '0');
		std::cout << "\n1a 1" << (x.size() == 3 ? "M" : x + "K") << std::endl;
		runq1(10000 * pow(10, i-1), genome, stk);
	}
}

void q2(int c, std::string genome, std::vector<Data>& stk) {
	for (int i = 1; i <= 3; i++) {
		string x(i, '0');
		std::cout << "\n1b 1" << (x.size() == 3 ? "M" : x + "K") << std::endl;
		runq2(10000 * pow(10,i-1), genome, stk);
	}
}

int main(int argc, char* argv[]) {
	assert(argc >= 3);
#ifdef TEST
#else
	std::ifstream ifs(argv[1]);
	assert(ifs.is_open());
#endif
	std::string genome;
	for (std::string str; std::getline(ifs, str); ) {
		if (str[0] == '>') continue;
		genome += str;
	}
	
	std::vector<Data> stk;
	if (argc == 4) {
		if (strcmp(argv[3], "q1") == 0) {
			q1(1, genome, stk);
		}
		else if (strcmp(argv[3], "q2") == 0) {
			q2(1, genome, stk);
		} else if (strcmp(argv[3], "q3") == 0) {
			ProcessDataset(genome, argv[2], 1000);
			//ProcessDataset(genome, argv[2], 10000);
			//ProcessDataset(genome, argv[2], 100000);
		} else if (strcmp(argv[3], "q4") == 0) {
			std::cout << "q4\n";
			vector<vector<int>> s;
			vector<vector<std::string>> t;
			std::string s1 = "AGCGTATCGCATGCATTCGCGCATAAGCTAG",
						s2 = "TCTCTGGAGCGGGCTTCGTATATGCTAAAGC";

			auto p = Blast_DB::query(s1, s2, &s, &t);
			std::cout << "Score: " << p.first << '\n';
			
			auto p2 = p.second;
			std::cout << p2.first << '\n';
			for (int i = 0; i < p2.first.size(); i++) {
				if (p2.first[i] != '-' && p2.second[i] != '-') {
					if (p2.first[i] != p2.second[i])
						std::cout << "x";
					else
						std::cout << "|";
				} else {
					std::cout << " ";
				}
			}
			std::cout << '\n';
			std::cout << p2.second << "\n\n";
			std::cout << "Scoring array:\n";
			print(s);
			std::cout << "Traceback array:\n";
			print(t);
			
		}
	} else {
		std::cout << "Invalid number of arguments.\n";
	}
}
