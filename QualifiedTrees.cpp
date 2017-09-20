#include "QualifiedTrees.h"
#include "main.h"
//#include <ctime>
//#include <cstdio>


QualifiedTreeList::QualifiedTreeList()
{
    m_listQualifiedTrees.clear();
    m_ticker = 0;
    initializeHistStatis();
    m_mapLayerWidthMatrix.clear();
    m_achieveValueTreeSpace = false;
    m_specificTreeNum = 0;
    m_foundIsomorphic = false;
    m_firstTree.delRoot();
    m_valueTreeSpace = 0;
}

QualifiedTreeList::~QualifiedTreeList()
{
    std::list<Tree>::iterator  iter = m_listQualifiedTrees.begin();
    for (; iter!= m_listQualifiedTrees.end(); iter++)
    {
        delete iter->m_root;
        iter->m_root = NULL;
    }
    m_listQualifiedTrees.clear();

    if (!m_firstTree.isEmptyTree()) {
        m_firstTree.delRoot();
    }
}

bool QualifiedTreeList::compareTree(const Tree& first, const Tree& second)
{
    if (first.m_property.m_nodesNum != second.m_property.m_nodesNum) return first.m_property.m_nodesNum > second.m_property.m_nodesNum;
    else if (first.m_property.m_layers != second.m_property.m_layers) return first.m_property.m_layers < second.m_property.m_layers;
    else if (first.m_property.m_maxEpsilon != second.m_property.m_maxEpsilon) return first.m_property.m_maxEpsilon < second.m_property.m_maxEpsilon;
    else return (first.m_property.m_filename.compare(second.m_property.m_filename) < 0);
}

void QualifiedTreeList::printTreeStatistic()
{

    if (0 == m_listQualifiedTrees.size()) return;

    std::string & originalVAFFilename = csvFile.m_inputFilename;

    //get time string
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char str[160];
    int wLen = sprintf(str, "_Stat_%d%02d%02d_%02d%02d%02d.csv",
                       timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,
                       timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
    str[wLen] = '\0';

    int strLen = originalVAFFilename.length();
    std::string filename = originalVAFFilename.substr(0,strLen-4)+str;


    FILE* pFile = fopen(filename.c_str(), "w");

    //print table head
    std::string csvHead("OriginalVAFFile,RootGene,TreeFilename,"
                        "CandidatesGeneNum,TreeGeneNum,inputTriangleEps(%),"
                        "minEps(%),aveEps(%),maxEps(%),totalEps(%),"
                        "Layers,Richness,criteriaValue\n");
    strLen = (int)csvHead.size();
    fwrite(csvHead.c_str(), sizeof(char), strLen, pFile);

    for (std::list<Tree>::iterator iter = m_listQualifiedTrees.begin();
            iter != m_listQualifiedTrees.end(); iter++)
    {

        std::string outLine = iter->m_property.m_originalVAFFilename+ ","
                              + csvFile.getVAFNode(iter->m_property.m_rootID).m_gene+ ","
                              +iter->m_property.m_filename+ ",";

        strLen = sprintf(str, "%d,%d,%d,%f,%.2f,%f,%d,%d,%.4f,%.4f\n",
                         iter->m_property.m_originalVAFNodesNum,iter->m_property.m_nodesNum,iter->m_property.m_inputEpsilon,
                         iter->m_property.m_minEpsilon,iter->m_property.m_averageEpsilon,iter->m_property.m_maxEpsilon,
                         iter->m_property.m_totalEpsilon,
                         iter->m_property.m_layers, iter->m_property.m_richness,iter->m_property.m_criteriaValue);
        str[strLen] = '\0';
        outLine.append(std::string(str));
        fwrite(outLine.c_str(), sizeof(char), outLine.length(), pFile);
    }

    fflush(pFile);
    fclose(pFile);

    printf("Output: %s\n", filename.c_str());

}

int QualifiedTreeList::getTicker()
{
    m_ticker++;
    if (m_ticker >= 10001) m_ticker = 0;
    return m_ticker;
}

void QualifiedTreeList::checkAddTree(Tree* pTree)
{
    pTree->computerTreeProperty(getTicker());

    if (constraint.m_printAllTrees != 0 && 0 == constraint.m_numQualifiedTrees){
        if ( (constraint.m_delete123 == 0 && 1 == pTree->m_property.m_richness) || constraint.m_delete123 > 0){
                       pTree->print(pTree->m_property.m_filename);
                       if (0 == m_specificTreeNum) {
                           m_firstTree = *pTree;
                           pTree->m_root = NULL;
                        }
                        else{
                           isIsomorphicWithFirstTree(pTree);
                        }
                       ++m_specificTreeNum;
        }

    }
    if (0 != constraint.m_analyzeStatistic) {
        updateLayerWidthMatrix(pTree->getLayerWidthVector());
    }

    //statistics all trees
    updateHistoricalStatis(pTree);

    if (0 == constraint.m_numQualifiedTrees){
        pTree->delRoot();
    }
    else {
        if ((int)m_listQualifiedTrees.size() >= constraint.m_numQualifiedTrees  && !compareTree(*pTree,m_listQualifiedTrees.back()))
        {
            pTree->delRoot();
        }
        else {
            if (!IsomorphicInList(pTree) && 1 == pTree->m_property.m_richness )
            {
                if ( (int)m_listQualifiedTrees.size() < constraint.m_numQualifiedTrees)
                {
                    m_listQualifiedTrees.push_back(*pTree);
                    pTree->m_root = NULL;
                    m_listQualifiedTrees.sort(compareTree);

                }
                else
                {
                    if (compareTree(*pTree,m_listQualifiedTrees.back()))
                    {
                        m_listQualifiedTrees.back().delRoot();
                        m_listQualifiedTrees.pop_back();
                        m_listQualifiedTrees.push_back(*pTree);
                        pTree->m_root = NULL;
                        m_listQualifiedTrees.sort(compareTree);
                    }
                    else //directly discard this tree;
                    {
                        pTree->delRoot();
                    }

                }

                doesAchieveValueTreeSpace();

            }
            else{
                updateHistIsomorphy(pTree);
                pTree->delRoot();

            }


        }

    }



    //print

    if (constraint.isElapsed() || doesAchieveValueTreeSpace())
    {
        if ((0 == constraint.m_analyzeStatistic || 2 == constraint.m_analyzeStatistic) && (constraint.m_numQualifiedTrees > 0) )
        {
            std::list<Tree>::iterator iter = m_listQualifiedTrees.begin();
            while (iter != m_listQualifiedTrees.end())
            {
                if (! iter->m_property.m_printed)
                {
                    iter->print(iter->m_property.m_filename);
                    iter->m_property.m_printed = true;
                }
                iter++;
            }
        }

        //print statistics table;
        printTreeStatistic();

        printHistoricalStatis();

        if (0 != constraint.m_analyzeStatistic) printLayerWidthMatrix();

        constraint.resetBeginningTime();

    }

}

bool QualifiedTreeList::IsomorphicInList(Tree* pTree)
{
    std::list<Tree>::iterator iter = m_listQualifiedTrees.begin();
    while (iter != m_listQualifiedTrees.end())
    {
        if (iter->isIsomorphic(pTree))
        {
            return true;
        }
        iter++;
    }
    return false;
}


void QualifiedTreeList::initializeHistStatis(){
    HistStatis temp;
    temp.m_numTree = 0;
    //temp.m_proportionTree = 0.0;
    temp.m_totalLayers = 0;
    temp.m_numIsomorphic = 0;
   // temp.m_proportionIsomorphic = 0;
    for (int i=0;i<11;i++) m_vectorHistStatis.push_back(temp);

}

int QualifiedTreeList::getPosInVectorFromRichness(float const richness){
    return int (richness*10);
}

 void QualifiedTreeList::updateHistoricalStatis(Tree* pTree){
    int index = getPosInVectorFromRichness(pTree->m_property.m_richness);
    m_vectorHistStatis.at(index).m_numTree++;
    m_vectorHistStatis.at(index).m_totalLayers += pTree->m_property.m_layers;
}

 void QualifiedTreeList::updateHistIsomorphy(Tree* pTree){
    int index = getPosInVectorFromRichness(pTree->m_property.m_richness);
    m_vectorHistStatis.at(index).m_numIsomorphic++;
 }


 void QualifiedTreeList::printHistoricalStatis(){

    //if (0 == m_listQualifiedTrees.size()) return;

    std::string& originalVAFFilename = csvFile.m_inputFilename;

    //get time string
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char str[40];
    int wLen = sprintf(str, "_Hist_%d%02d%02d_%02d%02d%02d.csv",
                       timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,
                       timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
    str[wLen] = '\0';

    int strLen = originalVAFFilename.length();
    std::string filename = originalVAFFilename.substr(0,strLen-4)+str;


    FILE* pFile = fopen(filename.c_str(), "w");

    //print table head
    std::string csvHead("Richness_Rank,NumTrees,TreeProportion(%),AverageLayers,NumIsomorphic,IsomorphicProportion(%)\n");
    strLen = (int)csvHead.size();
    fwrite(csvHead.c_str(), sizeof(char), strLen, pFile);

    int vSize = m_vectorHistStatis.size();

    //get total trees
    unsigned long nTotalTrees = 0;
    unsigned long nTotalIsomorphicTrees = 0;
    for (int i=0; i<vSize;i++){
        nTotalTrees += m_vectorHistStatis.at(i).m_numTree;
        nTotalIsomorphicTrees += m_vectorHistStatis.at(i).m_numIsomorphic;
    }

    char rowText[120];
    for (int i=0; i<vSize;i++)
    {

        float averLayers = 0;
        if (0 != m_vectorHistStatis.at(i).m_numTree)
        {
            averLayers = (float)m_vectorHistStatis.at(i).m_totalLayers*1.0/m_vectorHistStatis.at(i).m_numTree;

        }

        if (10 != i){
            strLen = sprintf(rowText, "%.1f__%.3f,%lu,%.1f,%.2f,%lu,%.1f\n",
                             i*0.1,
                             (i+1)*0.1-0.001,
                             m_vectorHistStatis.at(i).m_numTree,
                             100*(m_vectorHistStatis.at(i).m_numTree*1.0/nTotalTrees),
                             averLayers,
                             m_vectorHistStatis.at(i).m_numIsomorphic,
                             100*(m_vectorHistStatis.at(i).m_numIsomorphic*1.0/nTotalTrees) );

        }
        else
        {

            strLen = sprintf(rowText, "1.00,%lu,%.1f,%.2f,%lu,%.1f\n",
                             m_vectorHistStatis.at(i).m_numTree,
                             100*(m_vectorHistStatis.at(i).m_numTree*1.0/nTotalTrees),
                             averLayers,
                             m_vectorHistStatis.at(i).m_numIsomorphic,     100*(m_vectorHistStatis.at(i).m_numIsomorphic*1.0/nTotalTrees));
        }
        rowText[strLen] = '\0';
        fwrite(rowText, sizeof(char), strLen, pFile);
    }

    //print sum rows;
    strLen = sprintf(rowText,"Sum,%lu,100,,%lu,%.2f\n",
                     nTotalTrees,nTotalIsomorphicTrees,100*(nTotalIsomorphicTrees*1.0/nTotalTrees));
    rowText[strLen] = '\0';
    fwrite(rowText, sizeof(char), strLen, pFile);

    fflush(pFile);
    fclose(pFile);

    printf("Output: %s\n", filename.c_str());


 }


void QualifiedTreeList::updateLayerWidthMatrix(std::vector<int> const& treeLayerWidthVector){
    int size = int(treeLayerWidthVector.size());

    for(int i = 0; i<size; i++){
        m_mapLayerWidthMatrix[i][treeLayerWidthVector.at(i)] +=1;
    }
}

void QualifiedTreeList::printLayerWidthMatrix(){

    if (0 == m_mapLayerWidthMatrix.size()) return;

    std::string &originalVAFFilename = csvFile.m_inputFilename;

    //get time string
    time_t rawtime;
    struct tm * timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    char str[80];
    int wLen = sprintf(str, "_LayerWidthMatrix_%d%02d%02d_%02d%02d%02d.csv",
                       timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,
                       timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
    str[wLen] = '\0';

    int strLen = originalVAFFilename.length();
    std::string filename = originalVAFFilename.substr(0,strLen-4)+str;

    FILE* pFile = fopen(filename.c_str(), "w");

    //get maximum width
    int maxWidth = 0;
    std::map<int,std::map<int,long>>::iterator rowIter = m_mapLayerWidthMatrix.begin();
    while (rowIter != m_mapLayerWidthMatrix.end()){
        std::map<int,long>::iterator colEndIter =  --(rowIter->second.end());
        int tempWidth = colEndIter->first;
        if (tempWidth > maxWidth) maxWidth = tempWidth;
        rowIter++;
    }

    //print table head and define columnVector
    std::vector<int> columnVector;
    std::string csvHead("Width,1,");
    columnVector.push_back(1);
    for (int i=2;i<=maxWidth;i+=2) {
        char a[4];
        if (i != maxWidth) wLen = sprintf(a,"%d,",i);
        else wLen = sprintf(a,"%d\n",i);
        columnVector.push_back(i);
        a[wLen] = '\0';
        csvHead = csvHead+a;
    }
    strLen = (int)csvHead.size();
    fwrite(csvHead.c_str(), sizeof(char), strLen, pFile);

    //print row
    int colSize = columnVector.size();
    int nRow = m_mapLayerWidthMatrix.size();
    for(int i=0;i<nRow;i++){

        std::string rowText;
        rowText.clear();

        char a[16];
        wLen = sprintf(a,"Generation_%d,",i+1);
        a[wLen] = '\0';
        rowText.append(a);

        for(int j= 0; j<colSize; j++){
            if (0 == m_mapLayerWidthMatrix[i].count(columnVector.at(j))){
                if (j != colSize -1) rowText.append("0,");
                else rowText.append("0\n");
            }
            else{
                if (j != colSize -1) wLen = sprintf(a,"%lu,",m_mapLayerWidthMatrix[i][columnVector.at(j)]);
                else wLen = sprintf(a,"%lu\n",m_mapLayerWidthMatrix[i][columnVector.at(j)]);
                a[wLen] = '\0';
                rowText.append(a);
            }
        }

        fwrite(rowText.c_str(), sizeof(char), rowText.length(), pFile);

    }

    fflush(pFile);
    fclose(pFile);

    printf("Output: %s\n", filename.c_str());


}

bool QualifiedTreeList::doesAchieveValueTreeSpace(){
    m_achieveValueTreeSpace = (0 == constraint.m_numQualifiedTrees)? false : (m_listQualifiedTrees.size() >= constraint.m_numQualifiedTrees);
    return m_achieveValueTreeSpace;

}
bool QualifiedTreeList::doesAchieveSpecificTreeNum(){
    if (constraint.m_printAllTrees <= 0 ) return false;
    else {
        return m_specificTreeNum >= constraint.m_printAllTrees;
    }

}

bool QualifiedTreeList::isIsomorphicWithFirstTree(Tree* pTree){
      bool result = m_firstTree.isIsomorphic(pTree);
      if (result){
          printf("Information: %s and %s are isomorphic trees.\n",
                 m_firstTree.m_property.m_filename.c_str(), pTree->m_property.m_filename.c_str());
          m_foundIsomorphic = true;
      }
      return result ;

}

bool QualifiedTreeList::exitRandomConstruction(){
     return m_foundIsomorphic && m_specificTreeNum >= m_valueTreeSpace && m_specificTreeNum >= constraint.m_printAllTrees;
}


