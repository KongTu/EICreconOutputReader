#!/bin/bash


root -b -q src/pfRICH_QA.cxx+\(\"./fileLists/fileList_${1}x${2}.list\",\"DIS_${1}x${2}\",${1},${2}\)


