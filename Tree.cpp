#include "Tree.h"
//#include <list>
#include <cmath>
//#include <cstdio>
#include "main.h"
//#include <string>
#include <iostream>

TreeNode::TreeNode()
{
    m_nodeID = -1;
    m_nodeAlt = 0;
    m_layerIndex = 0;
    m_leftChild = NULL;
    m_rightChild = NULL;
    m_treeID.clear();
}

//nested delete all its children
TreeNode::~TreeNode()
{
    if (NULL != m_leftChild)
    {
        delete m_leftChild;
        m_leftChild = NULL;
    }

    if (NULL != m_rightChild)
    {
        delete m_rightChild;
        m_rightChild = NULL;
    }

}

TreeNode& TreeNode::operator = (TreeNode const& rhs)
{
    if (this != &rhs)
    {
        this->m_nodeID = rhs.m_nodeID;
        this->m_nodeAlt = rhs.m_nodeAlt;
        this->m_layerIndex = rhs.m_layerIndex;
        if (NULL != rhs.m_leftChild)
        {
            this->m_leftChild = new TreeNode();
            *(this->m_leftChild) = *(rhs.m_leftChild);
        }
        if (NULL != rhs.m_rightChild)
        {
            this->m_rightChild = new TreeNode();
            *(this->m_rightChild) = *(rhs.m_rightChild);
        }
    }
    return *this;
}


void  TreeNode::constructTreeID(std::string myTreeID){
    this->m_treeID = myTreeID;
    if (NULL != this->m_leftChild && NULL != this->m_rightChild)
    {
        m_leftChild->constructTreeID(myTreeID+"0");
        m_rightChild->constructTreeID(myTreeID+"1");
    }
}

//if 2 comparing nodes in 2 trees have same alt vaf instead of NodeID, it is isomorphic
bool TreeNode::isIsomorphicNode(TreeNode const* rhs){
    if (rhs->m_nodeAlt != this->m_nodeAlt){
        return false;
    }
    else{
        if (NULL != m_leftChild && NULL != m_rightChild
            && NULL != rhs->m_leftChild && NULL != rhs->m_rightChild)
        {
            return (m_leftChild->isIsomorphicNode(rhs->m_leftChild) && m_rightChild->isIsomorphicNode(rhs->m_rightChild))
            || (m_leftChild->isIsomorphicNode(rhs->m_rightChild) && m_rightChild->isIsomorphicNode(rhs->m_leftChild));
        }
        else return (NULL == m_leftChild && NULL == m_rightChild
                 && NULL == rhs->m_leftChild && NULL == rhs->m_rightChild);//all children are NULL


   }

}

Tree::Tree()
{
    m_root = NULL;

}

Tree::~Tree()
{
    if (NULL != m_root)
    {
        delete m_root;
        m_root = NULL;
    }
}

void Tree::setRoot(TreeNode* root)
{
    m_root = root;
}

void Tree::delRoot()
{
    if (NULL != m_root)
    {
        delete m_root;
        m_root = NULL;
    }
    delProperty();
}

bool Tree::isEmptyTree()
{
    return NULL == m_root;

}

void Tree::constructTreeID(){
    m_root->constructTreeID(std::string("1"));
}



// BFS search tree nodes
TreeNode* Tree::getNode(int const nodeID)
{
    if (NULL == m_root) return NULL;

    TreeNode* result = NULL;
    std::list<TreeNode*>  heap;
    heap.push_back(m_root);
    while (0 != heap.size())
    {
        TreeNode* cur = heap.front();
        heap.pop_front();
        if (nodeID ==cur->m_nodeID)
        {
            result = cur;
            break;
        }
        else
        {
            if (NULL != cur->m_leftChild) heap.push_back(cur->m_leftChild);
            if (NULL != cur->m_rightChild) heap.push_back(cur->m_rightChild);
        }
    }
    heap.clear();
    return result;
}


void Tree::addChildren(TreeNode* parent,TreeNode* leftChild, TreeNode* rightChild)
{
    parent->m_leftChild = leftChild;
    parent->m_rightChild = rightChild;

}

Tree& Tree::operator = (Tree const& rhs)
{
    if (this != &rhs)
    {
        this->m_root = rhs.m_root;
        this->m_property = rhs.m_property;
    }
    return *this;
}

// return AverageEps, parameter return Total Eps;

 void Tree::print(const std::string& filename){
    //constructTreeID();

    FILE* pFile = fopen(filename.c_str(), "w");

	//print table head
	int strLen = 0;
	std::string csvHead("TreeID,Layer,GENE,CHROM,POS,REF,ALT,DEPTH,REF_READ_NUM,ALT_READ_NUM,REF_FREQ,ALT_FREQ\n");
	strLen = (int)csvHead.size();
	fwrite(csvHead.c_str(), sizeof(char), strLen,pFile);

	//print data using BFS search
	std::string outLine;
	std::list<TreeNode*> heap;
	heap.push_back(m_root);//push back root index
	while (0 != heap.size()) {
		TreeNode* node = heap.front();
		heap.pop_front();
		VAF vaf = csvFile.getVAFNode(node->m_nodeID);
		if (NULL != node->m_leftChild) heap.push_back(node->m_leftChild);
		if (NULL != node->m_rightChild) heap.push_back(node->m_rightChild);

		//construct line text;
		outLine.clear();
		char nucleotides[4];
		nucleotides[3] = '\0';
		sprintf(nucleotides, "%c,%c", vaf.m_reference, vaf.m_alternative);
		outLine = node->m_treeID + ","
			+ std::to_string(node->m_layerIndex) + ","
			+ vaf.m_gene + ","
			+ vaf.m_chromosome + ","
			+ vaf.m_position + ","
			+ std::string(nucleotides) + ","
			+ std::to_string(vaf.m_depth) + ","
			+ std::to_string(vaf.m_refRead) + ","
			+ std::to_string(vaf.m_altRead) + ","
			+ std::to_string(vaf.m_refFreq*1.0/100) + ","
			+ std::to_string(vaf.m_altFreq*1.0/100) + "\n";
		strLen = (int)outLine.size();
		fwrite(outLine.c_str(), sizeof(char), strLen, pFile);
	}
	fflush(pFile);
	fclose(pFile);

 }

 void Tree::computerTreeProperty(int ticker)
{
    m_property.m_originalVAFFilename = csvFile.m_inputFilename;
    m_property.m_originalVAFNodesNum = (int)csvFile.m_vafVector.size();
    m_property.m_inputEpsilon = constraint.m_triConstraintEps;
    m_property.m_rootID = m_root->m_nodeID;

    m_property.m_richness = 0;
    m_property.m_totalEpsilon = 0;
    m_property.m_averageEpsilon = 0.0;
    m_property.m_minEpsilon = 100;
    m_property.m_maxEpsilon = 0.0;
    m_property.m_layers = 0;
    m_property.m_nodesNum = 0;
    m_property.m_criteriaValue = 1000;
    int triagNum = 0;

    //get Tree Epsilon
    if (NULL != m_root){
       std::list<TreeNode*> heap;
       heap.push_back(m_root);
       while (0 != heap.size()){
            TreeNode* node = heap.front();
            heap.pop_front();
            m_property.m_nodesNum++;

            if (node->m_layerIndex > m_property.m_layers) m_property.m_layers = node->m_layerIndex;

            if (NULL != node->m_leftChild && NULL != node->m_rightChild){
                triagNum++;
                int tempEps = abs(node->m_leftChild->m_nodeAlt + node->m_rightChild->m_nodeAlt - node->m_nodeAlt);

                if (tempEps < m_property.m_minEpsilon) m_property.m_minEpsilon = tempEps;
                if (tempEps > m_property.m_maxEpsilon) m_property.m_maxEpsilon = tempEps;
                m_property.m_totalEpsilon += tempEps;

                heap.push_back(node->m_leftChild);
                heap.push_back(node->m_rightChild);
            }

       }
       if (0 != triagNum) m_property.m_averageEpsilon = m_property.m_totalEpsilon*1.0/triagNum;
    }

    m_property.m_richness = m_property.m_nodesNum*1.0/m_property.m_originalVAFNodesNum;
    m_property.m_criteriaValue = (m_property.m_totalEpsilon*1.0/100+1.0)/(m_property.m_richness*m_property.m_richness);

    //get time and construct filename
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char str[40];
    int wLen = sprintf(str, "_T_E%d%%_%d%02d%02d_%02d%02d%02d_%d.csv",
                  m_property.m_inputEpsilon,
                  timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,
                  timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec,
                  ticker);
    str[wLen] = '\0';

    int length = m_property.m_originalVAFFilename.length();
    m_property.m_filename = m_property.m_originalVAFFilename.substr(0,length-4)+str;

    //m_property.m_numIsomorph = 0;
    m_property.m_printed = false;

    constructTreeID();
    if (1 == constraint.m_analyzeStructure){
        analyzeTreeStructure();
        printTreeStructuresStatistics();
    }

}


void Tree::analyzeTreeStructure(){
     if (NULL != m_root){
       //initialize tree structure vector
       LayerStructure temp;
       temp.m_balance = 0;
       temp.m_numLeaf = 0;
       temp.m_width = 0;
       m_property.m_layerStructVector.clear();
       for (int i=0;i<m_property.m_layers;i++){
            temp.m_layerIndex = i+1;
            m_property.m_layerStructVector.push_back(temp);
       }

       //BFS to get layer-width and numLeaf
       std::list<TreeNode*> heap;
       heap.push_back(m_root);
       while (0 != heap.size()){
            TreeNode* node = heap.front();
            heap.pop_front();
            m_property.m_layerStructVector.at(node->m_layerIndex-1).m_width++;

            if (NULL != node->m_leftChild && NULL != node->m_rightChild){
                heap.push_back(node->m_leftChild);
                heap.push_back(node->m_rightChild);
            }
            else{ //Leaf node
                m_property.m_layerStructVector.at(node->m_layerIndex-1).m_numLeaf++;
            }

       }

       //compute the balance and treeWidth
       m_property.m_treeWidth = 0;
       m_property.m_layerStructVector.at(0).m_balance = 1;
       for (int i=1;i<m_property.m_layers;i++){
            LayerStructure & preLayer =  m_property.m_layerStructVector.at(i-1);
            LayerStructure & curLayer =  m_property.m_layerStructVector.at(i);
            curLayer.m_balance = curLayer.m_width*1.0/(2*preLayer.m_width);

            if (curLayer.m_width > m_property.m_treeWidth) m_property.m_treeWidth = curLayer.m_width;
       }
    }

}

//this function must use after constructTreeID()
std::vector<int> Tree::getLayerWidthVector(){
    int totalLayers = m_property.m_layers;
    std::vector<int> vectorLayerWidth;
    vectorLayerWidth.clear();
    for(int i=0;i<totalLayers; i++) vectorLayerWidth.push_back(0);


    std::list<TreeNode*> heap;
    heap.push_back(m_root);
    while (0 != heap.size())
    {
        TreeNode* node = heap.front();
        heap.pop_front();
        int layer = int(node->m_treeID.length())-1;
        vectorLayerWidth.at(layer) +=1;

        if (NULL != node->m_leftChild && NULL != node->m_rightChild)
        {
            heap.push_back(node->m_leftChild);
            heap.push_back(node->m_rightChild);
        }
    }
    return vectorLayerWidth;
}

void Tree::delProperty(){
    m_property.m_filename.clear();
    m_property.m_richness = 0;
    m_property.m_originalVAFFilename.clear();
    m_property.m_originalVAFNodesNum = 0;
    m_property.m_rootID = -1;
    m_property.m_inputEpsilon = 0;
    m_property.m_totalEpsilon = 0;
    m_property.m_averageEpsilon = 0;
    m_property.m_minEpsilon = 0;
    m_property.m_maxEpsilon = 0;
    m_property.m_layers = 0;
    m_property.m_nodesNum = 0;
    m_property.m_criteriaValue = 0;

    //m_property.m_numIsomorph = 0;
    m_property.m_printed = false;
}

bool Tree::isIsomorphic(const Tree* other){
    if (m_property.m_layers != other->m_property.m_layers
        || m_property.m_nodesNum != other->m_property.m_nodesNum
        || m_property.m_rootID != other->m_property.m_rootID)
    {
        return false;
    }
    else {
        return m_root->isIsomorphicNode(other->m_root);

    }


}

void Tree::printTreeStructuresStatistics(){

    //Only print richness = 1's tree
    if (m_property.m_originalVAFNodesNum != m_property.m_nodesNum) return;

    //get filename
    int strLen = m_property.m_originalVAFFilename.length();
    std::string filename = m_property.m_originalVAFFilename.substr(0,strLen-4)+"_GenerationStructure.csv";

    //try read file, if it does not exist, create it.
    FILE* pFile = fopen(filename.c_str(), "r");
    if (NULL == pFile)
    {
        //create a table head file
        pFile = fopen(filename.c_str(), "w");
        std::string csvHead("TreeFilename,TreeDepth,TreeWidth,_Generation,Width,Skewness,Balance,...\n");
        strLen = (int)csvHead.size();
        fwrite(csvHead.c_str(), sizeof(char), strLen, pFile);
        printf("Continue updating Output: %s\n", filename.c_str());
    }
    fclose(pFile);


    //formally append file
    pFile = fopen(filename.c_str(), "a");

    std::string outText = m_property.m_filename;

    char tempText[20];
    strLen = sprintf(tempText,",%d,%d,",m_property.m_layers,m_property.m_treeWidth);
    tempText[strLen] = '\0';
    outText.append(tempText);

    //print layer structure
    int size = (int)m_property.m_layerStructVector.size();
    for(int i=0; i<size; i++)
    {
        LayerStructure& layStruct = m_property.m_layerStructVector.at(i);
        if (i == size -1){
             strLen = sprintf(tempText,"_%d,%d,%d,%.4f\n",
                            layStruct.m_layerIndex,layStruct.m_width,layStruct.m_numLeaf,layStruct.m_balance);

        }
        else{
           strLen = sprintf(tempText,"_%d,%d,%d,%.4f,",
                            layStruct.m_layerIndex,layStruct.m_width,layStruct.m_numLeaf,layStruct.m_balance);

        }
        tempText[strLen] = '\0';
        outText.append(tempText);

    }
    strLen = outText.length();
    fwrite(outText.c_str(), sizeof(char), strLen, pFile);

    fflush(pFile);
    fclose(pFile);
}





