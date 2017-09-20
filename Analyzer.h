#ifndef ANALYZER_H
#define ANALYZER_H
#include "CSVFile.h"
#include "Constraints.h"
#include <map>
#include <list>
#include "Tree.h"

struct NodePair{
    int m_node1;
    int m_alt1;
    int m_node2;
    int m_alt2;
};

struct TriangleCandidates {
    int m_altParent;

    int m_fragmentStatus;
    //  0:  not a fragment;
    // 10:  left child has been used in a path;
    // 01:  right child has been used in a path;
    // 11:  both left and right children have been used in a path;
    //  2:  fragment Node has been used in tree construction process;

    std::list<NodePair> m_nodePairList;
};

class Analyzer
{
public:
    Analyzer();
    ~Analyzer();
    Tree m_tree;


    void setRelatedClassPara(CSVFile *pCsvFile,Constraints* pConstraints);
    void constructNodeMap();
    void printTrianglesStatistics();
    void printSameVAFGeneNameList();
    void printTDS1();

    void constructRandomTree();

protected:


private:
    CSVFile* m_pCsvFile;
    Constraints* m_pConstraints;
    int m_ticker;
    std::map<int, TriangleCandidates> m_nodeMap;


    void eraseParentChildrenInMap(std::map<int, TriangleCandidates>* pNodeMap,
                                 int const parentID,
                                 int const leftChildID,
                                 int const rightChildID);
    void updateTDS1fromGeneFragment();
    void updateTDS1fromGenePath();

    //check whether Waiting Expand Nodes is Greater than unused gene fragment
    bool areNecessaryFragmentsUsed(std::list<int>& heap, std::map<int, TriangleCandidates> & nodeMap);
    double computeGeneTreeSpace();
    double computeValueTreeSpace();
    int getNumValuePair(std::list<NodePair>& nodePairList);
    int findChildIDInPairList(const int childFreq, TriangleCandidates& triangleCandidates);
    int findChildIDInUnspecificParentGroup(const int parentFreq, const int childFreq);
    int matchChildIDinSpecificParentGroup(const int parentID, const int childFreq);
    void deleteNochildIDinSelfParentGroup(const int childID, std::list<NodePair>& nodePairList);
    void deleteTaggedChildIDinOtherParentGroup(const int parentID, const int childID);
    void printPathVector(std::vector<int> const & pathVector);



};

#endif // ANALYZER_H
