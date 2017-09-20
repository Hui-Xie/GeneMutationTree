cmake_minimum_required(VERSION 3.5)
project(MuTree)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Analyzer.cpp
    Analyzer.h
    Constraints.cpp
    Constraints.h
    CSVFile.cpp
    CSVFile.h
    main.cpp
    main.h
    QualifiedTrees.cpp
    QualifiedTrees.h
    Tree.cpp
    Tree.h
    VAF_CSVTree_struct.cpp
    VAF_CSVTree_struct.h)

add_executable(MuTree ${SOURCE_FILES})