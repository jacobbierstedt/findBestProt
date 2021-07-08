#include <iostream>
#include "orfFind.hpp"



int main(int argc, const char *argv[]) {

	//A short nucleotide sequence
	std::string seq = "ATGCCCTCCGGAAAAAATTCATTTTAAGCCTGCTAA";

	//The reverse complement of the same sequence
	std::string rcseq = "TTAGCAGGCTTAAAATGAATTTTTTCCGGAGGGCAT";

	Orf o(seq);
	Orf rco(seq);

	std::cout << "Predicted protein sequence from seq:\t" << o.chooseLongestProtein() << std::endl;
	std::cout << "Predicted protein sequence from rcseq:\t" << rco.chooseLongestProtein() << std::endl;

	return 0;
}
