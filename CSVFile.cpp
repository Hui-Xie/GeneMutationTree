#include "CSVFile.h"
//#include "VAF_CSVTree_struct.h"
//#include <string>
//#include <cstdio>
#include <cstring>
//#include <cstdlib>
#include <iostream>
#include <cmath>
//#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>


#define LineLength 300



using namespace std;


bool fragmentCompare(const Fragment & first, const Fragment & second){
	if (first.m_parent != second.m_parent) return first.m_parent > second.m_parent;
	else return first.m_child >= second.m_child;
}
bool pathCompare(const std::vector<int>& first, const std::vector<int>& second){

    std::vector<int>::const_iterator iter1 = first.begin();
    std::vector<int>::const_iterator iter2 = second.begin();
    while(iter1 != first.end() && iter2 != second.end()){
        if (*iter1 != *iter2) return   *iter1 > *iter2;
        else {
            ++iter1;
            ++iter2;
        }

    }
    if (iter1 == first.end() && iter2 == second.end()){
        //same path
        return true;
    }
    else if (iter1 != first.end()){
        return true;
    }
    else{
        return false;
    }

}

CSVFile::CSVFile()
{
	m_pInputFile = NULL;
	m_pOutputFile = NULL;
	m_fileType = 0;
	m_fragmentFile = "";
	m_fragmentVector.clear();
	m_pathMatrix.clear();
}


CSVFile::~CSVFile()
{
	if (NULL != m_pInputFile) {
        fclose(m_pInputFile);
        m_pInputFile =NULL;
	}
	if (NULL != m_pOutputFile) {
        fclose(m_pOutputFile);
        m_pOutputFile = NULL;
	}
}



int CSVFile::determineFileType(){
	m_fileType = 0;
	char firstLine[LineLength];
	int temp = readFirstLine(firstLine, LineLength);
	if (-1 == temp) return m_fileType;


	int nCommas = countCommas(firstLine);
	if (VAF9 == nCommas && 0 == getFirstItem(firstLine).compare("GENE")) m_fileType = VAF9;
	else if (VAF10 == nCommas && 0 == getFirstItem(firstLine).compare("GENE")) m_fileType = VAF10;
	else if (TREE11 == nCommas && 0 == getFirstItem(firstLine).compare("TreeID")) m_fileType = TREE11;
	else m_fileType = 0;

	return m_fileType;
}

int CSVFile::countCommas(char* lineText){
	char *p = lineText;
	int counter = 0;
	while ('\0' != *p){
		if (',' == *p)  counter++;
		p++;
	}
	return counter;
}

std::string CSVFile::getFirstItem(char* lineText){
	char *p = lineText;
	int counter = 0;
	while (',' != *p && '\0'!=*p){
		counter++;
		p++;
	}
	counter= counter+1;
	char* str = new char[counter];
	for (int i = 0; i < counter - 1; i++){
		*(str + i) = *(lineText + i);
	}
	*(str + counter - 1) = '\0';
	std::string item(str);
	delete[] str;
	//item = item.TrimLeft().TrimRight();
	return item;
}

void CSVFile::removeAllSpaces(char* buffer, int nSize){
	char* iter = buffer;
	int nCount = 0;
	while ('\0' == *iter ){
		if (' ' == *iter){
			memcpy(iter, iter + 1, nSize - nCount - 1);
		}
		nCount++;
		iter++;
	}
}

//result: p+result will point to The Carriage
int  CSVFile::getOffsetToCarriage(char* p){
	char* iter = p;
	int nCount = 0;
	while ('\n' != *iter && '\0' != *iter){
		nCount++;
		iter++;
	}
	return nCount;
}

void  CSVFile::initializeVAFNode(VAF& node){

	node.m_reference = ' ';
	node.m_alternative = ' ';
	node.m_depth = 0;
	node.m_refRead = 0;
	node.m_altRead = 0;
	node.m_refFreq = 0;
	node.m_altFreq = 0;
	//node.m_2pq = 0;

}

void CSVFile::initializeBinNode(Node& node){
	node.m_index = -1;
	node.m_layer = -1;
	node.m_treeID = "";
	node.m_parent = -1;
	node.m_lchild = -1;
	node.m_rchild = -1;
}

//return the pointer in the next line
char* CSVFile::parseVAFLine(char* p, VAF& node){

	if ('\0' == *p) return p;

	initializeVAFNode(node);

	char* iter = p;
	int nItem = 0;
	while (true){
		int itemOffset = 0;
		char *itemBeginning = iter;
		while (',' != *iter && '\n' != *iter && '\0' != *iter){
			itemOffset++;
			iter++;
		}
		nItem++;
		int itemLength = itemOffset + 1;
		char* item = new char[itemLength];
		memcpy(item, itemBeginning, itemOffset);
		item[itemLength - 1] = '\0';

		switch (nItem){
		case 1:
			node.m_gene = std::string(item);
			break;
		case 2:
			node.m_chromosome = std::string(item);
			break;
		case 3:
			node.m_position = std::string(item);
			break;
		case 4:
			node.m_reference = *item;
			break;
		case 5:
			node.m_alternative = *item;
			break;
		case 6:
			node.m_depth = atoi(item);
			break;
		case 7:
			node.m_refRead = atoi(item);
			break;
		case 8:
			node.m_altRead = atoi(item);
			break;
		case 9:
			//node.m_refFreq = roundf((float)atof(item), 2);
			node.m_refFreq = round((float)atof(item)*100);
			break;
		case 10:
			//node.m_altFreq = roundf((float)atof(item), 2);
			node.m_altFreq = round((float)atof(item)*100);
			break;
		case 11:
			//node.m_2pq = roundf((float)atof(item), 2);
			break;
		default:
			cout<<"There is an error in Parsing VAF line. item: "<<item<<endl;

		}
		delete[] item;

		if ('\n' == *iter)  return ++iter; //point to next line
		else if ('\0' == *iter)  return iter;
		else iter++; //current iter pointing to comma
	}

}


char* CSVFile::parseTreeLine(char* p, VAF& vafNode, Node& binNode){

	if ('\0' == *p) return p;

	initializeVAFNode(vafNode);
	initializeBinNode(binNode);

	char* iter = p;
	int nItem = 0;
	while (true){
		int itemOffset = 0;
		char *itemBeginning = iter;
		while (',' != *iter && '\n' != *iter && '\0' != *iter){
			itemOffset++;
			iter++;
		}
		nItem++;
		int itemLength = itemOffset + 1;
		char* item = new char[itemLength];
		memcpy(item, itemBeginning, itemOffset);
		item[itemLength - 1] = '\0';

		switch (nItem){
		case 1:
			binNode.m_treeID = std::string(item);
			break;
		case 2:
			binNode.m_layer = atoi(item);
			break;
		case 3:
			vafNode.m_gene = std::string(item);
			break;
		case 4:
			vafNode.m_chromosome = std::string(item);
			break;
		case 5:
			vafNode.m_position = std::string(item);
			break;
		case 6:
			vafNode.m_reference = *item;
			break;
		case 7:
			vafNode.m_alternative = *item;
			break;
		case 8:
			vafNode.m_depth = atoi(item);
			break;
		case 9:
			vafNode.m_refRead = atoi(item);
			break;
		case 10:
			vafNode.m_altRead = atoi(item);
			break;
		case 11:
			//vafNode.m_refFreq = roundf((float)atof(item),2);
			vafNode.m_refFreq = round((float)atof(item)*100);
			break;
		case 12:
			//vafNode.m_altFreq = roundf((float)atof(item),2);
			//binNode.m_altFreq = roundf((float)atof(item),2);
			vafNode.m_altFreq = round((float)atof(item)*100);
			binNode.m_altFreq = round((float)atof(item)*100);
			break;
		default:
			cout<<"There is an error in Parsing Tree line. item:"<<item<<endl;
		}
		delete[] item;

		if ('\n' == *iter)  return ++iter; //point to next line
		else if ('\0' == *iter)  return iter;
		else iter++; //current iter pointing to comma
	}
}

//return the effective length of first line
int CSVFile::readFirstLine(char* pLine,const int size){
	FILE* pFile = NULL;
	pFile = fopen(m_inputFilename.c_str(), "rb");
	if (NULL == pFile) {
		cout<<"Sorry. Can not open file. If it is open, please close it and retry."<<endl;
		return -1;
	}
	fread(pLine,sizeof(char),size,pFile);
	fclose(pFile);
	//pFile = NULL;


	int length = 0;

	bool isFirstLine = true;
	for (int i = 0; i < size; i++){
		if ('\n' == *(pLine + i) && isFirstLine) {
			isFirstLine = false;
			length = i + 1;
			continue;
		}
		if (!isFirstLine)  *(pLine + i) = '\0';
	}
	return length;
}

float CSVFile::roundf(float x, int nDecimal)
{
	int ploid = 1;
	for (int i = 0;i < nDecimal;i++) {
		ploid = ploid * 10;
	}
	return round(x*ploid)/ploid;
}

int CSVFile::readFileData()
{
	int fileType = determineFileType();
	if (VAF9 != fileType && VAF10 != fileType){
        cout<<"Error: File type is not correct. Please check it."<<endl;
        return -1;
	}
	m_pInputFile = fopen(m_inputFilename.c_str(), "rb");
	if (NULL == m_pInputFile) {
		cout<<"Sorry. Can not open VAF file. If it is open, please close it and retry."<<endl;
		return -11;
	}

	fseek(m_pInputFile,0,SEEK_END);
    int nFileSize = ftell(m_pInputFile);
	fseek(m_pInputFile,0,SEEK_SET);

	char* buffer = new char[nFileSize+1];
	int readSize = fread(buffer, sizeof(char), nFileSize, m_pInputFile);
	fclose(m_pInputFile);
	m_pInputFile = NULL;

    *(buffer + readSize) = '\0';
    char* p = buffer;
    removeAllSpaces(buffer, nFileSize + 1);
    int lineOffset = 0;
    lineOffset = getOffsetToCarriage(p);
    p += lineOffset + 1; //let p point to first Non-head data row;


    while ('\0' != *p)
    {
        VAF newNode;
        p = parseVAFLine(p, newNode);
        m_vafVector.push_back(newNode);
    }

    deleterZeroNode();

    m_vafVector.shrink_to_fit();
    sortVAFVector();

    delete[] buffer;
    return 0;
}

void CSVFile::addFragment(const int parent, const int child){
	Fragment temp;
	temp.m_parent = parent;
	temp.m_child = child;
	m_fragmentVector.push_back(temp);
}


int CSVFile::readFragmentFile(){
	std::ifstream fragmentFileStream;
	std::filebuf * fb = fragmentFileStream.rdbuf();

	if (NULL == fb->open (m_fragmentFile.c_str(),std::ios_base::in)){
		cout<<"Sorry. Can not open fragment file. If it is open, please close it and retry."<<endl;
		return -11;
	}

	std::string line;
	std::string parent;
	std::string child;
	int parentFreq =0;
	int childFreq =0;
	int nLine = 0;  //record number of line
	int nPath = 0;
    int nChildren = 0;
	while(std::getline(fragmentFileStream,line))
	{
		++nLine;
		while (' ' == line.at(0)) line.erase(0,1);
		if ('#'==line.at(0)) continue;
		std::stringstream  lineStream(line);
		std::getline(lineStream,parent,',');
		parentFreq = stoi(parent);
		++nPath;
        nChildren  = 0;
		//get fragment
        while(std::getline(lineStream,child,','))
		{
			try{
				childFreq = stoi(child);
			}
			catch(const std::invalid_argument& ia){
				printf("\n****Error: Maybe the fragment.csv file has redundancy comma at line end.\n");
				printf("    suggest: exit program to check fragment.csv file again.\n");
	     	}

			if (childFreq*2 >= parentFreq) 	addFragment(parentFreq,childFreq);
			else {
				printf("\n****Error: in line #%d, Parent = %d, child = %d. It should be child*2 >= parent.\n", nLine,parentFreq,childFreq);
				printf("       suggest: exit program to check fragment.csv file again.\n");
			}
			parentFreq = childFreq;
            ++nChildren;
     	}


	}
	fb->close();
	std::sort(m_fragmentVector.begin(),m_fragmentVector.end(), fragmentCompare);
    printf("Progress: program gets %d gene fragments from %d gene paths.\n", m_fragmentVector.size(),nPath);
	return 0;
}

int CSVFile::readPathFile(){
	std::ifstream pathFileStream;
	std::filebuf *fb = pathFileStream.rdbuf();

	if (NULL == fb->open(m_fragmentFile.c_str(), std::ios_base::in)) {
		cout << "Sorry. Can not open fragment file. If it is open, please close it and retry." << endl;
		return -11;
	}

	std::string line;
	std::string child;
	int childFreq = 0;
	int nLine = 0;  //record number of line
	int nPath = 0;
	while (std::getline(pathFileStream, line)) {
		++nLine;
		if (line.size() < 3){
			printf("Error: at row %d of path File file,incorrect input format: %s\n", nLine, line.c_str());
			continue;
		}
		while (' ' == line.at(0)) line.erase(0, 1);
		int temp =line.size()-1;
		if ('#' == line.at(0)
			|| !(line.at(0)>= '0' && line.at(0) <= '9'
				 && line.at(temp)>= '0' && line.at(temp) <= '9') ) continue;

		//get path vector
		std::stringstream pathStream(line);
		std::vector<int> pathVector;
		pathVector.clear();
		while (std::getline(pathStream, child, ',')) {
			try {
				childFreq = stoi(child);
			}
			catch (const std::invalid_argument &ia) {
				printf("\n****Error: Maybe the path file has redundancy comma at line %d end.\n",nLine);
				printf("    suggest: exit program to check fragment.csv file again.\n");
			}
			pathVector.push_back(childFreq);
		}
		m_pathMatrix.push_back(pathVector);
		++nPath;
	}
	fb->close();
	std::sort(m_pathMatrix.begin(), m_pathMatrix.end(), pathCompare);
	printf("Progress: program gets %d gene paths.\n", nPath);
	return 0;

}

void CSVFile::deleterZeroNode(){
    int nOriginalSize = (int)m_vafVector.size();
    std::vector<VAF>::iterator iter = m_vafVector.begin();
    while (iter != m_vafVector.end()){
        if (0 == iter->m_altFreq) iter = m_vafVector.erase(iter);
        else iter++;
    }
    int nNewSize = (int)m_vafVector.size();
    printf("Progress: Deleted %d zero nodes, remains %d nodes.\n",nNewSize-nOriginalSize,nNewSize);
}


//in descending order of m_altFreq
void CSVFile::sortVAFVector(){
	int size = (int) m_vafVector.size();
	for (int i = 0; i < size - 1; i++){
		int max = i;
		for (int j = i + 1; j < size; j++){
			if (m_vafVector.at(j).m_altFreq > m_vafVector.at(max).m_altFreq) max = j;
		}

		if (max != i){
			VAF temp;
			temp = m_vafVector.at(i);
			m_vafVector.at(i) = m_vafVector.at(max);
			m_vafVector.at(max) = temp;
		}
	}

	//build the m_mapIDGenename
	m_mapIDGenename.clear();
	for (int i=0;i<size;i++) m_mapIDGenename[i] = m_vafVector.at(i).m_chromosome;
}

VAF CSVFile::getVAFNode(int const geneID){
    return m_vafVector.at(geneID);
}


