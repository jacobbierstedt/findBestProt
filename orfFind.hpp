#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>



class Orf {
private:
	std::string _seq;
	std::string orfFind(std::string nuc);
	std::string getLongestOrf(std::string prot);
	std::string revComp(std::string nuc);
	std::string codonToAminoAcid(std::string codon);

public:
	Orf(std::string &seq) : _seq(seq) {};
	virtual ~Orf() {};
	void printSeq();
	std::string chooseLongestProtein();
};



void Orf::printSeq() {
std::cout << _seq << std::endl;
std::cout << revComp(_seq) << std::endl;
}



std::unordered_map<std::string, std::string> genCode({
																			{"ATG" , "M"},
																			{"TTT" , "F"},
																			{"TTC" , "F"},
																			{"TTA" , "L"},
																			{"TTG" , "L"},
																			{"TCT" , "S"},
																			{"TCC" , "S"},
																			{"TCA" , "S"},
																			{"TCG" , "S"},
																			{"TAT" , "Y"},
																			{"TAC" , "Y"},
																			{"TAA" , "*"},
																			{"TAG" , "*"},
																			{"TGT" , "C"},
																			{"TGC" , "C"},
																			{"TGA" , "*"},
																			{"TGG" , "W"},
																			{"CTT" , "L"},
																			{"CTC" , "L"},
																			{"CTA" , "L"},
																			{"CTG" , "L"},
																			{"CCT" , "P"},
																			{"CCC" , "P"},
																			{"CCA" , "P"},
																			{"CCG" , "P"},
																			{"CAT" , "H"},
																			{"CAC" , "H"},
																			{"CAA" , "Q"},
																			{"CAG" , "Q"},
																			{"CGT" , "R"},
																			{"CGC" , "R"},
																			{"CGA" , "R"},
																			{"CGG" , "R"},
																			{"ATT" , "I"},
																			{"ATC" , "I"},
																			{"ATA" , "I"},
																			{"ACT" , "T"},
																			{"ACC" , "T"},
																			{"ACA" , "T"},
																			{"ACG" , "T"},
																			{"AAT" , "N"},
																			{"AAC" , "N"},
																			{"AAA" , "K"},
																			{"AAG" , "K"},
																			{"AGT" , "S"},
																			{"AGC" , "S"},
																			{"AGA" , "R"},
																			{"AGG" , "R"},
																			{"GTT" , "V"},
																			{"GTC" , "V"},
																			{"GTA" , "V"},
																			{"GTG" , "V"},
																			{"GCT" , "A"},
																			{"GCC" , "A"},
																			{"GCA" , "A"},
																			{"GCG" , "A"},
																			{"GAT" , "D"},
																			{"GAC" , "D"},
																			{"GAA" , "E"},
																			{"GAG" , "E"},
																			{"GGT" , "G"},
																			{"GGC" , "G"},
																			{"GGA" , "G"},
																			{"GGG" , "G"}
																});


// take codon and return amino acid
std::string Orf::codonToAminoAcid(std::string codon){
		return genCode[codon];
}


// Take nucleotide seq and return translation and longest orf
std::string Orf::orfFind(std::string nuc){
	int nlen = nuc.size() - (nuc.size() % 3);
	std::string prot("");
	for(int i=0, inc = 3; i < nlen ; i += inc){
		prot = prot+codonToAminoAcid(nuc.substr(i,3));
	}
	return getLongestOrf(prot);
}


// returns index of longest string in vector of strings
size_t maxOrfLength(std::vector<std::string> orfs){
	size_t max = 0;
	size_t maxIndex = 0;
	for(size_t i = 0; i < orfs.size(); i++){
		if(orfs[i].size() > max){
			max = orfs[i].size();
			maxIndex = i;
		}
	}
	return maxIndex;
}


// Takes protein sequence string and returns longest ORF
std::string Orf::getLongestOrf(std::string prot){
	std::vector<std::string> orfv;
	for(int i = 0; i < prot.size(); i++){
		std::string pep;
		if(prot[i] == 'M'){
			pep = prot.substr(i);
			orfv.push_back( pep.substr(0,pep.find('*')+1) );
		}
	}
	std::string longProt;
	if(!orfv.empty()){
		longProt = orfv[maxOrfLength(orfv)];
	} else {longProt = "";}
	return longProt;
}


// reverse complements nucleotide sequence
std::string Orf::revComp(std::string nuc){
	std::string rc("");
	for(size_t i=0;i<nuc.size();i++){
		if(nuc[i] == 'A'){
			rc = rc + 'T';
		}else if(nuc[i] == 'T'){
			rc = rc + 'A';
		}else if(nuc[i] == 'G'){
			rc = rc + 'C';
		}else if(nuc[i] == 'C'){
			rc = rc + 'G';
		}else{
			rc = rc + 'N';
		}
	}
	std::reverse(rc.begin(), rc.end());
	return rc;
}


// translate nucleotide sequence to 6 frames, choose longest orf beginning with M
// and ending with a stop codon.
std::string Orf::chooseLongestProtein(){
	std::vector<std::string> frames = {
		_seq,
		_seq.substr(1),
		_seq.substr(2),
		revComp(_seq),
		revComp(_seq).substr(1),
		revComp(_seq).substr(2)
	};

	std::vector<std::string> bestORFs;
	for(size_t i=0;i<frames.size();i++){
		bestORFs.push_back(orfFind(frames[i]));
	}
	std::string bestPep = bestORFs[maxOrfLength(bestORFs)];
	return bestPep;
}
