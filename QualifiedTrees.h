#ifndef QUALIFIEDTREES_H
#define QUALIFIEDTREES_H
#include "Tree.h"
#include <list>
#include <vector>
#include <map>


//historical Statistics of all computing trees
struct HistStatis{
    unsigned long m_numTree;
    //float m_proportionTree;
    unsigned long m_totalLayers;
    unsigned long m_numIsomorphic;
    //float m_proportionIsomorphic;
};


class QualifiedTreeList
{

    public:
        QualifiedTreeList();
        ~QualifiedTreeList();
        void checkAddTree(Tree* pTree);
        bool IsomorphicInList(Tree* pTree);
        bool doesAchieveValueTreeSpace();
        bool m_achieveValueTreeSpace;
        bool doesAchieveSpecificTreeNum();
        double m_valueTreeSpace;
        bool exitRandomConstruction();


    protected:

    private:
        std::list<Tree> m_listQualifiedTrees;
        int m_specificTreeNum;

        std::vector<HistStatis> m_vectorHistStatis;
        int getPosInVectorFromRichness(float const richness);
        void updateHistoricalStatis(Tree* pTree);
        void updateHistIsomorphy(Tree* pTree);
        void printHistoricalStatis();
        void initializeHistStatis();

        static bool compareTree(const Tree& first, const Tree& second);

        void printTreeStatistic();
        void printLayerWidthMatrix();

        int m_ticker;
        int getTicker();

        std::map<int,std::map<int,long>> m_mapLayerWidthMatrix;
        void updateLayerWidthMatrix(std::vector<int> const& treeLayerWidthVector);
        Tree m_firstTree;

        bool isIsomorphicWithFirstTree(Tree* pTree);
        bool m_foundIsomorphic;


};

#endif // QUALIFIEDTREES_H
