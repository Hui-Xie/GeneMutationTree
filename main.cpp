#include <iostream>
#include "main.h"
#include "Analyzer.h"


using namespace std;


//string curDir;

CSVFile csvFile;
Constraints constraint;
QualifiedTreeList qualifiedTrees;


void printParameterHints(){

    cout <<"=========================================================================================================\n"
    "====================================Mutation Tree Build Program==========================================\n\n"
    "Version: June 15th, 2016\n"
    "Description: Build gene mutation Value trees, which have difference in alternative CCF value instead of gene name.\n"
    "             Reconstruct gene mutation Value trees from gene fragments or paths information file.\n"
    "             Use Ctrl-C to exit this program.\n\n"
    "             Brother program FindIsomorphy will help you find isomorphic trees in current directory.\n"
    "Usage: MuTree [Option Option_parameter] ...\n\n"
    "EXAMPLES:\n"
    "MuTree -F filename -E 0.02 -R 0.40 -H 2.0 -Q 200\n"
    "MuTree -F VAF-TCGA-A5-A0G9.csv -E 1 -R 0.40 -H 0.2 -Q 400\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -E 0 -R 0.30 -H 0.3 -Q 50 -A 1\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -E 0 -R 0.30 -H 8 -Q 0 -A 1 -S 1\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -E 0 -R 0.30 -H 2 -P -1 \n"
    "MuTree -F VAF-SIMU-NUM-55.csv -E 0 -R 0.30 -H 0.02 -P 0 -D 1\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -G GeneFragment.csv -E 0 -R 0.30 -H 3 -P -1 -D 1\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -G GenePath.csv -E 0 -R 0.30 -H 3 -P 10000 -D 0\n"
    "MuTree -F VAF-SIMU-NUM-55.csv -G GenePath.csv -E 0 -R 0.30 -H 3 -P -1 -D 0\n\n"

    "OPTIONS:\n"
    "-F VAF filename;\n"
    "-G gene fragment or path filename, which is used to reconstruct gene tree from gene fragment or path information;\n"
    "-D delete all possible VAF 1+1,2+1,1+2 combinations to reduce computation;\n"
    "   0: default, not delete; 1: delete;\n"
    "-E triangle epsilon percent in integer,e.g. 1 indicate 1%;\n"
    "-R richness of tree\n"
    "-H minimum printing interval(hour), e.g. 2 or 0.5\n"
    "-Q number of storing qualified trees\n"
    "-A analyze the statistics relations between generations and width of all constructing trees:\n"
    "   0:No analyze; 1:analyze without printing tree; 2: analyze with printing tree\n"
    "-S analyze tree's generation structure including its width, balance, and skewness\n"
    "   0:default no analysis;\n"
    "   1:analyze and print generation structure, and automatically set Q = 0\n"
    "-P print all or specific number random constructed trees, e.g. -1, 0, 100000\n"
    "   0: default, not print trees; -1: print all random trees; n: (n>0) print specific n trees;\n"
    "=========================================================================================================\n";

}

int main(int argc, char** argv)
{

    if (1 != argc%2 || argc < 3) {
         printParameterHints();
         return -1;
    }

    //Check input arguments
    for(int i = 1; i<argc-1; i= i+2 ){
        std::string argLeader(*(argv+i));

        if ( 0 == argLeader.compare("-F") || 0 == argLeader.compare("-f")){
            csvFile.m_inputFilename = std::string (*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-E") || 0 == argLeader.compare("-e")){
            constraint.m_triConstraintEps = int(atof(*(argv+i+1)));
        }
        else if ( 0 == argLeader.compare("-R") || 0 == argLeader.compare("-r")){
            constraint.m_leastRichness = (float)atof(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-H") || 0 == argLeader.compare("-h")){
            constraint.m_printInterval = (float)atof(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-Q") || 0 == argLeader.compare("-q")){
            constraint.m_numQualifiedTrees = atoi(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-A") || 0 == argLeader.compare("-a")){
            constraint.m_analyzeStatistic = atoi(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-S") || 0 == argLeader.compare("-s")){
            constraint.m_analyzeStructure = atoi(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-P") || 0 == argLeader.compare("-p")){
            constraint.m_printAllTrees = atoi(*(argv+i+1));
            //if (constraint.m_printAllTrees < 0){
            //    cout<<"=======Information: program will print huge number of trees in current dir.=============="<<std::endl;
            //}
         }
        else if ( 0 == argLeader.compare("-D") || 0 == argLeader.compare("-d")){
            constraint.m_delete123 = atoi(*(argv+i+1));
        }
        else if ( 0 == argLeader.compare("-G") || 0 == argLeader.compare("-g")){
            csvFile.m_fragmentFile = std::string (*(argv+i+1));
            constraint.m_reconstructFromFragment = 1;
        }
        else
            {
            printParameterHints();
            return -2;
        }
    }

    if (1 == constraint.m_analyzeStructure || constraint.m_printAllTrees != 0)
    {
        constraint.m_numQualifiedTrees = 0;
    }


    //read VAF file
    if (csvFile.m_inputFilename.empty()){
        cout<<"Error:Input Filename is empty."<<endl;
        printParameterHints();
        return -3;
    }
    if (0 != csvFile.readFileData()){
        cout<<"Error: read VAF file error. program exit."<< endl;
        return -4;
    }

    //read Gene fragment file
    if (constraint.m_reconstructFromFragment >=1 ){
        if (0 != csvFile.readPathFile()) return -11;
    }

    //Statistics Triangle
    Analyzer analyzer;
    analyzer.setRelatedClassPara(&csvFile, &constraint);
    analyzer.constructNodeMap();
    analyzer.printSameVAFGeneNameList();
    analyzer.printTrianglesStatistics();

    // user hits Ctrl-C to exit loop
    printf("\n=====If you want to stop this program, Please hit Ctrl-C to exit program=======\n");
    printf("=====If you do NOT want to stop this program, Please avoid to unintentionally hit Ctrl-C =========\n");
    while(1){
         analyzer.constructRandomTree();
         if (qualifiedTrees.doesAchieveSpecificTreeNum()) {
             printf("\n\n******=====Program achieves specific tree num. Exit.===========********\n");
             break;
         }
         if (qualifiedTrees.exitRandomConstruction()){
             printf("\n\n******====Program finds at least one isomorphic resulting tree with the first tree and achieve value tree space. Exit.=======\n");
             break;
         }
    }

    return 0;
}
