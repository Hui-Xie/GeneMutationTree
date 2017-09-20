#ifndef TREE_H
#define TREE_H
#include <string>
#include <list>
#include <vector>


class TreeNode
{
public:
    TreeNode();
    ~TreeNode();

    int m_nodeID;
    int m_nodeAlt;
    int m_layerIndex;
    std::string m_treeID; //something like "10101"
    TreeNode* m_leftChild;
    TreeNode* m_rightChild;
    TreeNode& operator = (TreeNode const& rhs);
    void constructTreeID(std::string myTreeID);
    bool isIsomorphicNode(TreeNode const* rhs);

};

struct LayerStructure
{
    int m_layerIndex;
    int m_width;
    float m_balance;
    int m_numLeaf;// the number of leave node(1%)
};

struct TreeProperty
{
    std::string m_filename;
    float m_richness;
    std::string m_originalVAFFilename;
    int m_originalVAFNodesNum;
    int m_rootID;
    int m_inputEpsilon;
    int m_totalEpsilon;
    float m_averageEpsilon;
    float m_minEpsilon;
    float m_maxEpsilon;
    int m_layers;
    int m_nodesNum;
    float m_criteriaValue;

    //int m_numIsomorph;
    bool m_printed;

    std::vector<LayerStructure> m_layerStructVector;
    int m_treeWidth;
};


class Tree
{
public:
    Tree();
    ~Tree();

    TreeNode* m_root;
    TreeProperty m_property;


    void addChildren(TreeNode* parent,TreeNode* leftChild, TreeNode* rightChild);
    void setRoot(TreeNode* root);
    void delRoot();
    bool isEmptyTree();
    TreeNode* getNode(int const nodeID);
    Tree& operator = (Tree const & rhs);
    void print(const std::string& filename);

    void delProperty();

    void computerTreeProperty(int ticker);
    bool isIsomorphic(const Tree* other);

    std::vector<int> getLayerWidthVector();

protected:
private:

    void constructTreeID();
    void analyzeTreeStructure(); //analyze layer width, balance, skewness etc
    void printTreeStructuresStatistics();

};

#endif // TREE_H
