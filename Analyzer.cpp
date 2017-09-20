#include "Analyzer.h"
#include <cmath>
#include <assert.h>
//#include "stdio.h"
//#include <ctime>
//#include <cstdlib>
#include "main.h"

Analyzer::Analyzer()
{
    m_ticker = 0;
    srand (time(NULL));


}

Analyzer::~Analyzer()
{

}

void Analyzer::setRelatedClassPara(CSVFile *pCsvFile,Constraints* pConstraints){
    m_pCsvFile = pCsvFile;
    m_pConstraints = pConstraints;
}

void Analyzer::constructNodeMap(){
    printf("Progress: start to construct a node map.\n");

    int const& eps = m_pConstraints->m_triConstraintEps;
    std::vector<VAF> &  vafVector = m_pCsvFile->m_vafVector;
    m_nodeMap.clear();

    int size = vafVector.size();
    printf("Original vafVector size = %d\n", size);
    printf("Triangle Epsilon: %d%%\n", eps);

    //constructIDGenenameMap
    for(int i = 0; i<size-2; i++)
    {

        for(int j= i+1; j<size-1; j++)
        {
            for (int k= j+1; k<size; k++)
            {
                int diff = abs(vafVector.at(j).m_altFreq + vafVector.at(k).m_altFreq - vafVector.at(i).m_altFreq);
                if ((diff <= eps) || (0 == eps && diff < 1e-4))
                {
                    NodePair nodePair;
                    nodePair.m_node1 = j;
                    nodePair.m_node2 = k;
                    nodePair.m_alt1 = vafVector.at(j).m_altFreq;
                    nodePair.m_alt2 = vafVector.at(k).m_altFreq;

                    if ((constraint.m_delete123 > 0)  && ( nodePair.m_alt1 +  nodePair.m_alt2 <= 3)) continue;

                    if (0 == m_nodeMap.count(i)) {
                        m_nodeMap[i].m_altParent = vafVector.at(i).m_altFreq;
                        m_nodeMap[i].m_fragmentStatus = 0;
                     }
                    m_nodeMap[i].m_nodePairList.push_back(nodePair);
                }
            }
        }

    }



    printf("============================================================\n");
    printf("Before considering gene fragments, possible Gene Tree Space:  %e \n", computeGeneTreeSpace());
    printf("Before considering gene fragments, possible Value Tree Space:  %e \n", computeValueTreeSpace());
    printf("============================================================\n");

    if (m_pConstraints->m_reconstructFromFragment >0 ) {
        //updateTDS1fromGeneFragment();
        updateTDS1fromGenePath();
        printf("============================================================\n");
        printf("After considering gene fragments, possible Gene Tree Space:  %e \n", computeGeneTreeSpace());
        qualifiedTrees.m_valueTreeSpace = computeValueTreeSpace();
        printf("After considering gene fragments, possible Value Tree Space (Incorrect reference):  %e \n",qualifiedTrees.m_valueTreeSpace);
        printf("Possible value tree space is INCORRECT because same edge may be put in different layers\n");
        printf("============================================================\n");

    }

    printf("Total parent nodes number: %d\n", (int)m_nodeMap.size());
    printTDS1();
    printf("Progress: Finished the TDS1 construction\n\n");
}


int Analyzer::getNumValuePair(std::list<NodePair>& nodePairList){
    std::list<NodePair>::iterator iter = nodePairList.begin();
    int num = 0;
    int alt1=0,alt2=0;
    while (iter != nodePairList.end()){
        if (0== num){
            alt1 = iter->m_alt1;
            alt2 = iter->m_alt2;
            ++num;
        }
        if (alt1 != iter->m_alt1 || alt2 != iter->m_alt2){
            alt1 = iter->m_alt1;
            alt2 = iter->m_alt2;
            ++num;
        }
        ++iter;
    }
    return num;

}


//calculate the possible genne tree space
double Analyzer::computeGeneTreeSpace(){
    double space = 1.0;
    std::map<int, TriangleCandidates>::iterator mapIter = m_nodeMap.begin();
    while (mapIter !=m_nodeMap.end()){
        space = space * (mapIter->second.m_nodePairList.size());
        ++mapIter;
    }
    return space;
}

//calculate the possible value space
double Analyzer::computeValueTreeSpace(){
    double space = 1.0;
    std::map<int, TriangleCandidates>::iterator mapIter = m_nodeMap.begin();
    while (mapIter !=m_nodeMap.end()){
        space = space * getNumValuePair(mapIter->second.m_nodePairList);
        ++mapIter;
    }
    return space;

}

void Analyzer::updateTDS1fromGeneFragment(){
    std::vector<Fragment>& fragmentVector = m_pCsvFile->m_fragmentVector;
    if (fragmentVector.size() == 0) return;
    std::map<int, TriangleCandidates>::iterator mapIter = m_nodeMap.begin();
    std::vector<Fragment>::iterator fragIter = fragmentVector.begin();
    while (mapIter != m_nodeMap.end() && fragIter != fragmentVector.end() ){
        if (mapIter->second.m_altParent != fragIter->m_parent){
            while (mapIter->second.m_altParent > fragIter->m_parent && mapIter != m_nodeMap.end()) ++mapIter;
            while (mapIter->second.m_altParent < fragIter->m_parent && fragIter != fragmentVector.end()) ++fragIter;
            if (mapIter == m_nodeMap.end() || fragIter == fragmentVector.end()) break;
        }
        mapIter->second.m_fragmentStatus = 1;
        std::list<NodePair>::iterator pairIter = mapIter->second.m_nodePairList.begin();
        bool findChild = false;
        while (pairIter != mapIter->second.m_nodePairList.end()){
            if (pairIter->m_alt1 == fragIter->m_child || pairIter->m_alt2 == fragIter->m_child ) {
                ++pairIter;
                findChild = true;
            }
            else {
                pairIter = mapIter->second.m_nodePairList.erase(pairIter);
            }
        }
        if (!findChild) {
            printf("*******Error: the parent-child (%d - %d) fragment does not exist in possible combinations.\n", fragIter->m_parent,fragIter->m_child);
            printf("       Suggest to exit program, and check the fragment.csv.\n");
        }
        ++mapIter;
        ++fragIter;
    }
 }


void Analyzer::updateTDS1fromGenePath(){
    std::vector< std::vector<int> > pathMatrix = m_pCsvFile->m_pathMatrix;
    std::vector< std::vector<int> >::iterator rowIter = pathMatrix.begin();
    while(rowIter != pathMatrix.end()){
        std::vector<int>& pathVector = *rowIter;
        int pathSize = pathVector.size();
        bool isAncestor = true;
        int childID = -1;
        for (int i=0; i<pathSize-1;++i){
             if (isAncestor) {
                 childID = findChildIDInUnspecificParentGroup(pathVector.at(i),pathVector.at(i+1));
                 if (-1 == childID ) {
                     printf("\n**** Error: Repeated or illegal edge %d-%d in the path: ", pathVector.at(i),pathVector.at(i+1));
                     printPathVector(pathVector);
                     printf("**** Suggest: check input fragment or path file to make assure there are no repeated physical edges and invalidate nodes.\n");
                     break;
                 }
                 else {
                     isAncestor = false;
                 }
             }
             else {
                 int parentID = childID;
                 childID  = matchChildIDinSpecificParentGroup(parentID,pathVector.at(i+1));
                 if (-1 == childID) {
                     isAncestor = true;
                     printf("\n**** Error: illegal node %d in the path: ", pathVector.at(i+1));
                     printPathVector(pathVector);
                     printf("**** Suggest: check input fragment or path file to make assure there are no repeated physical edges and invalidate nodes.\n");
                     break;

                 }
              }

        }
        ++rowIter;
    }
}

void Analyzer::printTrianglesStatistics(){
    //construct file name
    int length = m_pCsvFile->m_inputFilename.size();
    char epsStr[7];
    sprintf(epsStr,"%d",m_pConstraints->m_triConstraintEps);
    epsStr[6] = '\0';

    std::string filename = m_pCsvFile->m_inputFilename.substr(0,length-4)
                  +"_Triangles_E"+epsStr+"%.csv";


    FILE* pFile = fopen(filename.c_str(), "w");

	std::string tableHead("Gene,altFreq(%),Num_Triangles\n");
	fwrite(tableHead.c_str(), sizeof(char), tableHead.size(), pFile);
	std::string outLine;

	std::map<int, TriangleCandidates>::iterator mapIter = m_nodeMap.begin();
    while (mapIter != m_nodeMap.end()){

        outLine.clear();
        outLine.append(m_pCsvFile->getVAFNode(mapIter->first).m_gene);
        outLine.append(",");

        char str[50];
        sprintf(str, "%d,%d\n", mapIter->second.m_altParent, (int)mapIter->second.m_nodePairList.size());

        outLine.append(std::string(str));

        fwrite(outLine.c_str(), sizeof(char), outLine.size(), pFile);
        mapIter++;
    }

	fflush(pFile);
	fclose(pFile);
	//pFile = NULL;

	printf("Output: %s\n", filename.c_str());


}

 void Analyzer::printSameVAFGeneNameList(){
    //construct file name
    int length = m_pCsvFile->m_inputFilename.size();
    std::string filename = m_pCsvFile->m_inputFilename.substr(0,length-4)+"_SameVAFList.csv";


    FILE* pFile = fopen(filename.c_str(), "w");

	std::string tableHead("ALT_FREQ,Num,Gene......\n");
	fwrite(tableHead.c_str(), sizeof(char), tableHead.size(), pFile);

	std::vector<VAF>::iterator iter = m_pCsvFile->m_vafVector.begin();
	int preAltFreq = iter->m_altFreq;
	std::vector<std::string> sameVAFGeneVector;
	std::string outLine;
    while (iter != m_pCsvFile->m_vafVector.end()){


        if (iter->m_altFreq != preAltFreq){
            //print preAltFreq;
            outLine.clear();
            char str[15];
            int nSameNum = (int)sameVAFGeneVector.size();
            int writeLen = sprintf(str, "%d%%,%d,", preAltFreq,nSameNum);
            str[writeLen] = '\0';
            outLine.append(std::string(str));
            for (int i=0; i< nSameNum;i++) {
                if (i == nSameNum-1){
                    outLine.append(sameVAFGeneVector[i]+"\n");
                }
                else {
                    outLine.append(sameVAFGeneVector[i]+",");

                }

            }
            fwrite(outLine.c_str(), sizeof(char), outLine.size(), pFile);

            //Set new AltFreq;
            preAltFreq = iter->m_altFreq;
            sameVAFGeneVector.clear();
        }

        sameVAFGeneVector.push_back(iter->m_gene);

        iter++;
    }

    //print the last batch same VAF list
    outLine.clear();
    char str[15];
    int nSameNum = (int)sameVAFGeneVector.size();
    int writeLen = sprintf(str, "%d%%,%d,", preAltFreq,nSameNum);
    str[writeLen] = '\0';
    outLine.append(std::string(str));
    for (int i=0; i< nSameNum; i++)
    {
        if (i == nSameNum-1)
        {
            outLine.append(sameVAFGeneVector[i]+"\n");
        }
        else
        {
            outLine.append(sameVAFGeneVector[i]+",");

        }

    }
    fwrite(outLine.c_str(), sizeof(char), outLine.size(), pFile);
    sameVAFGeneVector.clear();

    fflush(pFile);
    fclose(pFile);
    //pFile = NULL;

    printf("Output: %s\n", filename.c_str());


 }

// Print Triangle Data Set 1
void Analyzer::printTDS1(){
    //construct file name
    int length = m_pCsvFile->m_inputFilename.size();
    std::string filename = m_pCsvFile->m_inputFilename.substr(0,length-4)+"_TDS1.csv";

    FILE* pFile = fopen(filename.c_str(), "w");

    std::string tableHead("Parent_ALT_FREQ(%),IsFragment,NumGenePairs,NumValuePairs,GenePair(%),......\n");
    fwrite(tableHead.c_str(), sizeof(char), tableHead.size(), pFile);

    std::map<int, TriangleCandidates>::iterator iter = m_nodeMap.begin();
    std::string outLine;
    char str[15];
    while (iter != m_nodeMap.end()){
        outLine.clear();
        int numGenePairs = iter->second.m_nodePairList.size();
        int numValuePairs = getNumValuePair(iter->second.m_nodePairList);
        int writeLen = sprintf(str, "%d,%d,%d,%d,", iter->second.m_altParent,iter->second.m_fragmentStatus,numGenePairs, numValuePairs);
        str[writeLen] = '\0';
        outLine.append(std::string(str));
        std::list<NodePair>::iterator listIter = iter->second.m_nodePairList.begin();
        while (listIter != iter->second.m_nodePairList.end()){
            writeLen = sprintf(str, "%d+%d", listIter->m_alt1,listIter->m_alt2);
            str[writeLen] = '\0';
            outLine.append(std::string(str));
            if (listIter != --(iter->second.m_nodePairList.end())){
                outLine.append(",");
            }
            ++listIter;
        }
        outLine.append("\n");
        fwrite(outLine.c_str(), sizeof(char), outLine.size(), pFile);

        ++iter;
    }
    fflush(pFile);
    fclose(pFile);
    printf("Output: %s\n", filename.c_str());

}

void Analyzer::eraseParentChildrenInMap(std::map<int, TriangleCandidates>* pNodeMap,
                                     int const parentID,
                                     int const leftChildID,
                                     int const rightChildID)
{
     std::map<int, TriangleCandidates>::iterator iterMap = pNodeMap->begin();
     for(;iterMap != pNodeMap->end();++iterMap){
          if (parentID == iterMap->first) {
                pNodeMap->erase(iterMap);
                break;
          }
     }

     iterMap = pNodeMap->begin();
     for(;iterMap != pNodeMap->end();++iterMap){
         std::list<NodePair>::iterator iterList = iterMap->second.m_nodePairList.begin();
         while(iterList != iterMap->second.m_nodePairList.end()){
             if(leftChildID == iterList->m_node1
                || rightChildID == iterList->m_node1
                || leftChildID == iterList->m_node2
                || rightChildID == iterList->m_node2 ){

                iterList = iterMap->second.m_nodePairList.erase(iterList);
             }
             else{
                ++iterList;

             }
         }

     }
 }


void Analyzer::constructRandomTree(){

    //copy node map;
    std::map<int, TriangleCandidates> nodeMap = m_nodeMap;

    //BFS construct tree;
    std::list<int> heap;

    //get maximum node name
    heap.push_back(0);

    //Clear tree for further constructing
    m_tree.delRoot();

    m_ticker++;
    if (m_ticker > 10000){
        m_ticker = 0;
        srand (time(NULL));
    }

    //BSF construct tree
    bool successfulConstruct = true;
    while(0 != heap.size()){

        if (m_pConstraints->m_reconstructFromFragment > 0 && !areNecessaryFragmentsUsed(heap, nodeMap)) {
            successfulConstruct = false;
            break;
        }


        int parent = heap.front();
        heap.pop_front();

        if (0 == nodeMap.count(parent)) continue;
        TriangleCandidates& triCandidates = nodeMap[parent];
        if (1 == triCandidates.m_fragmentStatus || 10 == triCandidates.m_fragmentStatus || 11 == triCandidates.m_fragmentStatus){
                  triCandidates.m_fragmentStatus = 2;
        }


        if (m_tree.isEmptyTree())
        {
            TreeNode* treeNode = new TreeNode();
            treeNode->m_nodeID = parent;
            treeNode->m_nodeAlt = triCandidates.m_altParent;
            treeNode->m_layerIndex = 1;
            m_tree.setRoot(treeNode);
        }

        int listSize = triCandidates.m_nodePairList.size();
        if (0 == listSize) continue;
        int randomPos = 0;
        if (listSize > 1)  randomPos = rand()% listSize;

        std::list<NodePair>::iterator iter= triCandidates.m_nodePairList.begin();
        for (int i=0;i<randomPos;++i){
            ++iter;
        }

        TreeNode* parentNode= m_tree.getNode(parent);

        TreeNode* leftChild = new TreeNode();
        leftChild->m_nodeID = iter->m_node1;
        leftChild->m_nodeAlt = iter->m_alt1;
        leftChild->m_layerIndex = parentNode->m_layerIndex +1;

        TreeNode* rightChild = new TreeNode();
        rightChild->m_nodeID = iter->m_node2;
        rightChild->m_nodeAlt = iter->m_alt2;
        rightChild->m_layerIndex = parentNode->m_layerIndex +1;

        m_tree.addChildren(parentNode, leftChild, rightChild);

        heap.push_back(leftChild->m_nodeID);
        heap.push_back(rightChild->m_nodeID);

        eraseParentChildrenInMap(&nodeMap,parent, leftChild->m_nodeID, rightChild->m_nodeID);
    }

    if (successfulConstruct){
        //check tree's property
        qualifiedTrees.checkAddTree(&m_tree);

    }
    else{
        m_tree.delRoot();
    }

}

bool Analyzer::areNecessaryFragmentsUsed(std::list<int>& heap, std::map<int, TriangleCandidates>& nodeMap){
    if (0 == heap.size()) return true;

    //get maximum VAF in heap;
    int maximumVAF = nodeMap[heap.front()].m_altParent;
    std::list<int>::iterator heapIter = heap.begin();
    ++heapIter;
    while (heapIter != heap.end() ){
        if (0 != nodeMap.count(*heapIter)) {
            int temp = nodeMap[*heapIter].m_altParent;
            if (temp > maximumVAF)  maximumVAF = temp;
         }
         ++heapIter;
    }

    //get maximum unused fragment in nodeMap
    int unusedFragmentParent = 0;
    std::map<int, TriangleCandidates>::iterator mapIter = nodeMap.begin();
    while (mapIter != nodeMap.end()){
        if (1 == mapIter->second.m_fragmentStatus  || 10 == mapIter->second.m_fragmentStatus  || 11 == mapIter->second.m_fragmentStatus ){
            unusedFragmentParent = mapIter->second.m_altParent;
            break;
        }
        ++mapIter;
    }

    return maximumVAF >= unusedFragmentParent;

}

int  Analyzer::findChildIDInPairList(const int childFreq, TriangleCandidates& triangleCandidates){
    std::list<NodePair>& nodePairList = triangleCandidates.m_nodePairList;
    const int parentFragmentStatus = triangleCandidates.m_fragmentStatus;
    int childID = -1;
    std::list<NodePair>::const_iterator iter = nodePairList.begin();
    while (iter != nodePairList.end()){
        if (childFreq == iter->m_alt1 && (1 == parentFragmentStatus || 0 == parentFragmentStatus) ) {
            childID =   iter->m_node1;
            if (1 == parentFragmentStatus) triangleCandidates.m_fragmentStatus = 11;
            if (0 == parentFragmentStatus) triangleCandidates.m_fragmentStatus = 10;
            break;
        }
        if (childFreq == iter->m_alt2 && (10 == parentFragmentStatus || 0 == parentFragmentStatus) ) {
            childID =   iter->m_node2;
            if (10 == parentFragmentStatus) triangleCandidates.m_fragmentStatus = 11;
            if (0 == parentFragmentStatus) triangleCandidates.m_fragmentStatus = 01;
            break;

        }
        ++iter;
    }
    return childID;

}

int  Analyzer::findChildIDInUnspecificParentGroup(const int parentFreq, const int childFreq){
    int childID = -1;
    int parentID = -1;
    // It must first seek 01 and 10 untagged node,then 00 untagged node
    std::map<int, TriangleCandidates>::iterator iter = m_nodeMap.begin();
    while (iter != m_nodeMap.end()){
        if (parentFreq == iter->second.m_altParent && (1 == iter->second.m_fragmentStatus || 10 == iter->second.m_fragmentStatus) )
        {
            childID = findChildIDInPairList(childFreq,iter->second);
            if (-1  != childID) {
                parentID = iter->first;
                break;
             }
         }
        ++iter;
    }

    if (-1 == childID){
        iter = m_nodeMap.begin();
        while (iter != m_nodeMap.end()){
            if (parentFreq == iter->second.m_altParent  && 0 == iter->second.m_fragmentStatus)
            {
                childID = findChildIDInPairList(childFreq,iter->second);
                if (-1  != childID) {
                    parentID = iter->first;
                    break;
                }
            }
            ++iter;
        }
    }

    // delete tagged nodeID in self and other parent groups
    if (-1 != childID){
        //delete nonChildID in current nodePairList
        deleteNochildIDinSelfParentGroup(childID,iter->second.m_nodePairList);

        //delete childID in non currentParent group;
        deleteTaggedChildIDinOtherParentGroup(parentID,childID);
    }
    return childID;
}

int Analyzer::matchChildIDinSpecificParentGroup(const int parentID, const int childFreq){
    int childID = -1;
    if ( 1 != m_nodeMap.count(parentID) )  return childID;

    std::list<NodePair>& nodePairList = m_nodeMap[parentID].m_nodePairList;
    int const parentFragmentStatus = m_nodeMap[parentID].m_fragmentStatus;
    std::list<NodePair>::iterator iter = nodePairList.begin();

    while(iter != nodePairList.end() ){
        if (childFreq == iter->m_alt1 && (01 == parentFragmentStatus || 0 ==  parentFragmentStatus) ) {
            childID =  iter->m_node1;
            if (1 == parentFragmentStatus)  m_nodeMap[parentID].m_fragmentStatus = 11;
            if (0 == parentFragmentStatus)  m_nodeMap[parentID].m_fragmentStatus = 10;
            break;
        }
        if (childFreq == iter->m_alt2 && (10 == parentFragmentStatus || 0 ==  parentFragmentStatus) ) {
            childID = iter->m_node2;
            if (10 == parentFragmentStatus)  m_nodeMap[parentID].m_fragmentStatus = 11;
            if (0 == parentFragmentStatus)  m_nodeMap[parentID].m_fragmentStatus = 01;
            break;
        }
        ++iter;

    }

    // delete tagged nodeID in self and other parent groups
    if (-1 != childID){
        //delete nonChildID in current nodePairList
        deleteNochildIDinSelfParentGroup(childID,nodePairList);

        //delete childID in non currentParent group;
        deleteTaggedChildIDinOtherParentGroup(parentID,childID);
    }
    else{
        printf("Can not find childID: parentFragmentStatus= %d, parentVAF = %d, childVAF=%d in matchChildIDinSpecificParentGroup\n",
               parentFragmentStatus, m_nodeMap[parentID].m_altParent,childFreq);
        printf("Suggest: please check the input fragment or path file to assure there is no repeated physical fragment.\n ");
    }

    return childID;

}


//delete nonChildID in current nodePairList
void Analyzer::deleteNochildIDinSelfParentGroup(const int childID, std::list<NodePair>& nodePairList){
    std::list<NodePair>::iterator iter = nodePairList.begin();
    while(iter != nodePairList.end()){
        if (childID != iter->m_node1 &&  childID != iter->m_node2) {
            if (nodePairList.size() >= 2) iter = nodePairList.erase(iter);
            else break;
        }
        else {
            ++iter;
        }
    }
}


void Analyzer::deleteTaggedChildIDinOtherParentGroup(const int parentID, const int childID){
    std::map<int, TriangleCandidates>::iterator mapIter = m_nodeMap.begin();
    while (mapIter != m_nodeMap.end()){
        if (parentID != mapIter->first){
            std::list<NodePair>& pairList = mapIter->second.m_nodePairList;
            std::list<NodePair>::iterator pairIter = pairList.begin();
            while (pairIter != pairList.end() ){
                if (childID == pairIter->m_node1 || childID == pairIter->m_node2){
                    if  (pairList.size() >= 2)  pairIter = pairList.erase(pairIter);
                    else break;
                }
                else {
                    ++pairIter;
                }
             }
        }
        ++mapIter;
    }

}

void Analyzer::printPathVector(std::vector<int> const & pathVector){
    std::vector<int>::const_iterator iter = pathVector.begin();
    while (iter != pathVector.end() ){
        if (iter != --(pathVector.end()))  printf("%d,",*iter);
        else printf("%d\n",*iter);
        ++iter;
    }
}




