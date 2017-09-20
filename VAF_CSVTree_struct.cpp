
#include "VAF_CSVTree_struct.h"

Node::Node()
{
	m_index = -1;
	m_layer = -1;
	m_treeID = "";
	m_parent = -1;
	m_lchild = -1;
	m_rchild = -1;
	m_altFreq = 0;
}

Node::~Node()
{
}

Node & Node::operator=(const Node & other)
{
	this->m_index = other.m_index;
	this->m_layer = other.m_layer;
	this->m_treeID = other.m_treeID;
	this->m_parent = other.m_parent;
	this->m_lchild = other.m_lchild;
	this->m_rchild = other.m_rchild;
	this->m_altFreq = other.m_altFreq;
	return *this;
}
