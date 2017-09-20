#ifndef _VAF_CSVTree_Struct_H_
#define _VAF_CSVTree_Struct_H_

#include <string>

struct VAF
{
	std::string m_gene;
	std::string m_chromosome;
	std::string m_position;
	char	m_reference;
	char	m_alternative;
	int 	m_depth;
	int 	m_refRead;
	int 	m_altRead;
	//float   m_refFreq;
	//float 	m_altFreq;
	int     m_refFreq; //original m_refFreq*100
	int     m_altFreq; //original m_altFreq*100
	//float  	m_2pq;    // m_refFreq* m_altFreq
};


class Node{
public:
	Node();
	~Node();

	int m_index;
	int m_layer;
	std::string m_treeID;
	int m_parent; //vector index in the  m_nodeVector, also the index in VAFVector
	int m_lchild; //vector index in the  m_nodeVector, also the index in VAFVector
	int m_rchild; //vector index in the  m_nodeVector, also the index in VAFVector
	int m_altFreq;
	Node& operator =(const Node& other);

};

#endif // _VAF_CSVTree_Struct_H_

