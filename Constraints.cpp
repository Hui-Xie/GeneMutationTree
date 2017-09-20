#include "Constraints.h"
//#include <string>



Constraints::Constraints()
{
	//m_rootConstraintEps = 0.02;
	m_triConstraintEps = 0.05;
	m_leastRichness = 0.60;
	m_printInterval = 24; //unit:hour
	time(&m_beginningTime);
	m_numQualifiedTrees = 300;
	m_type = 2;//default BFS
	m_analyzeStatistic = 0;
	m_analyzeStructure = 0;
	m_printAllTrees = 0;
	m_delete123 = 0;
	m_reconstructFromFragment = 0;

}


Constraints::~Constraints()
{
}

bool Constraints::isElapsed(){
    time_t now;
    time(&now);
    double tDiff = difftime(now,m_beginningTime);
    if (tDiff > m_printInterval*3600){
       time(&m_beginningTime);
       return true;
    }
    else{
       return false;
    }

}

void Constraints::resetBeginningTime(){
    time(&m_beginningTime);
}



//1 = DFS; 2 = BFS; 3=
std::string Constraints::getTypeName(){
	switch (m_type){
	case 1:
		return std::string("DFS");
	case 2:
		return std::string("BFS");
	default:
		return std::string("");
	}
	//return std::string("");

}
