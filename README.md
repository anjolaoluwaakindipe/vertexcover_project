# VERTEXCOVER_PROJECT
This is was a project undertaken for my ECE 650 class in order to understand and analyse different vertex cover algorithms. In order to, calculate the efficiency of all algorithms, each algorithm was placed in different threads and both computational time and approximation ratio was taken into consideration for each algorithm. Computational Time was calculated based on how long it took each algorithm to get a result and approximation ration being how correct the result was compared to the optimal minimum cover, thus the amount of nodes in the result divided by the optimal amount of nodes to cover the graph.
Different algorithms and their optimizations were used in the analysis and more could be read about them in the report.pdf file in this repository.

## INFORMATION BEFORE RUNNING THE PROJECT
- This project uses the pthread library to create multiple threads thus it is intended to be ran in a UNIX-like operating system. 
- This project is dependent on the minisat library so make sure to download it before running the project. The library can be gotten from this [here](https://github.com/agurfinkel/minisat)
- Graphs are intended to be inputed using the following structure:
    ```bash
        V 5
        E {<1,3>,<5,2>,<1,2>,<1,5>,<1,4>,<2,3>,<4,2>}
    ```
    where **"V"** is the number of vertexes and **"E"** are the edges. The Vertexes must range from 1 to n inclusive if the number of vertices is n.
- You can find examples from the test.txt file. This file can be piped directly into the program after building or you could write your own test cases for this project
- A csv file of the result of this project will be created under {BUILD_DIRECTORY}/example.csv. It will contain the following data of each graph passed into it:
1. Number of vertex
2. Type of algorithm used
3. Computational Time
4. Approximation Ratio (If it could not determine a ratio, it would give -1)

```csv 
5,CNF,220.783,-1
5,3CNF,285.407,-1
5,VC1,75.972,1
5,VC2,65.824,2
5,VC1-REF,50.293,1
5,VC2-REF,20.469,1
15,CNF,1.13759e+07,-1
15,3CNF,1.39465e+07,-1
15,VC1,46.988,1
15,VC2,31.77,1.5
15,VC1-REF,38.091,1
15,VC2-REF,53.551,1
```
- Note that the result of an algorithm will be passed to the standard output displaying the cover gotten from a graph containing the name of the algorithm and the nodes in considers a cover. E.g
```text
CNF-SAT-VC: 1,2
CNF-3-SAT-VC: 1,2
APPROX-VC-1: 1,2
APPROX-VC-2: 1,2,3,5
REFINED-APPROX-VC-1: 1,2
REFINED-APPROX-VC-2: 1,2
CNF-SAT-VC: 1,2,3,4,5
CNF-3-SAT-VC: 1,2,3,4,5
APPROX-VC-1: 1,2,3,4,5
APPROX-VC-2: 1,2,3,4,5,8,9,10
REFINED-APPROX-VC-1: 1,2,3,4,5
REFINED-APPROX-VC-2: 1,2,4,5,8,9
```

## BUILDING THE PROJECT
1. The project uses CMake as a build tool thus create a BUILD_DIRECTORY which could be called anything.
```bash
mkdir build
```
2. Navigate into that directory and inside the directory type "cmake" and the path where the CMakeLists.txt file is in.
```bash
cd build
cmake ../
```
3. After Cmake has finish generating the make files that will be used to build the project type **"make"**
```bash
make
```

## RUNNING THE PROJECT
1. While still in the BUILD_DIRECTORY you could either call the prjece650 and type each graph individually
```bash 
./prjece650
V 5
E {<1,3>,<5,2>,<1,2>,<1,5>,<1,4>,<2,3>,<4,2>}
```
2. You could also pipe the test.txt file in the project directory which has 100 graphs to generate graphs from (Note: Because two algorithms in this analysis has an exponential time increase in relation to the number of vertices expect the program to run longer for graphs with 20 vertices and more)
```bash
./prjece650 < ../test.txt
```
3. You can also create you own test case in a text file and pipe it through the program using method above.
```bash
./prjece650 < ../{yourtestcase}.txt
```

#### Note
Please note that the program always appends to the example.csv files, so remember to clear the result if you needed to.


## THINGS I LEARNT
1. Multi-threading in linux using the pthread library
2. Analysing efficiency of multiple vertex cover algorithms
3. Learning about SAT solvers and how it could be used to solve NP problems
4. Optimizing algorithms for correctness vs speed
5. Using CMake as build system
