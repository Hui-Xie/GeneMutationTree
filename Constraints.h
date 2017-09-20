#ifndef _Constraints_H_
#define _Constraints_H_
#include <string>
#include <ctime>



class Constraints
{
public:
	//float m_rootConstraintEps;
	int m_triConstraintEps;
	float m_leastRichness;
	float m_printInterval; //unit:hour
    bool isElapsed();//whether the time elapse beyond the interval
    int  m_numQualifiedTrees;
    void resetBeginningTime();
	int m_printAllTrees;
	int m_delete123;
	int m_reconstructFromFragment;

    //0:No analyze; 1:analyze without printing tree; 2: analyze with printing tree
    int  m_analyzeStatistic;

    //analyze tree's generation structure
    int  m_analyzeStructure;

	//1 = DFS; 2 = BFS; 3=
	int m_type;

public:
	Constraints();
	~Constraints();

	std::string getTypeName();
private:
    time_t m_beginningTime;

};

#endif

