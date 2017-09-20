
#ifndef _CSVFile_H_
#define _CSVFile_H_

#include "VAF_CSVTree_struct.h"

#include <vector>
#include <map>

#define VAF9  9
#define VAF10 10
#define TREE11 11


/*
return value: 0, Error file format

m_filetype: VAF9:VAF file without 2pq; 9 commas;
GENE,CHROM,POS,REF,ALT,DEPTH,REF_READ_NUM,ALT_READ_NUM,REF_FREQ,ALT_FREQ\n

m_filetype: VAF10:VAF file with 2pq; 10 commas;
GENE,CHROM,POS,REF,ALT,DEPTH,REF_READ_NUM,ALT_READ_NUM,REF_FREQ,ALT_FREQ,2pq\n

m_filetype: TREE11: Tree File; 11 commas;
TreeID,Layer,GeneName,Chromosome,Position,Reference,Alternative,Depth,RefRead,AltRead,RefFreq,AltFreq\n
*/

struct Fragment{
	int m_parent;
	int m_child;
};

bool fragmentCompare(const Fragment & first, const Fragment & second);

class CSVFile
{
public:
	std::string m_inputFilename;
	FILE* m_pInputFile;
	FILE* m_pOutputFile;
	//std::string m_outputFileName;
	std::string m_fragmentFile;

	int m_fileType;
	std::vector<VAF>  m_vafVector;
	std::map<int,std::string> m_mapIDGenename;
	std::vector<Fragment> m_fragmentVector;
	std::vector< std::vector<int> > m_pathMatrix;


public:
	CSVFile();
	~CSVFile();

	int readFileData();
	int determineFileType();
	void sortVAFVector();
	VAF getVAFNode(int const geneID);
	int readFragmentFile();
	int readPathFile();



private:
	int countCommas(char* lineText);
	std::string getFirstItem(char* lineText);
	int getOffsetToCarriage(char* p);
	void removeAllSpaces(char* buffer, int nSize);
	char* parseVAFLine(char* p, VAF& node);
	char* parseTreeLine(char* p, VAF& vafNode, Node& binNode);
	void initializeVAFNode(VAF& node);
	void initializeBinNode(Node& node);

	int readFirstLine(char* pLine, const int size);

	float roundf(float x, int nDecimal);
	void deleterZeroNode();

	void addFragment(const int parent, const int child);


};

#endif // _CSVFile_H_

