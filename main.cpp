// Compile with c++ a2ece650.cpp -std=c++11 -o a2ece650
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <regex>
#include <queue>
#include <set>
#include <limits>
#include <algorithm>
#include <minisat/core/Solver.h>

class GraphStructure
{
public:
    std::map<int, std::vector<int>> treeStructure;
    int numberOfVertices;
    std::vector<int> edgeList;
    std::map<int, int> treeNodeLevel;
    std::queue<int> nodeQueue;
    std::set<int> nodeSet;

    // methods
    void addToEdgeList(int node1, int node2);
    void generateGraphMap(int numberOfVertices);
    void matchNodeByEdges(std::string &args);
    std::vector<int> findPath(int currentNode, int startingNode, std::map<int, std::vector<int>> treeStructure, std::map<int, int> treeNodeLevel);
    void findPathBetweenDifferentNodes(int &fromNode, int &toNode, bool &found, std::vector<int> &resultPath);
    void findPathBetweenSameNodes(int &fromNode, int &toNode, bool &found, std::vector<int> &resultPath);
};

void GraphStructure::addToEdgeList(int node1, int node2)
{
    this->edgeList.push_back(node1);
    this->edgeList.push_back(node2);
}

void GraphStructure::generateGraphMap(int numberOfVertices)
{

    this->treeStructure.clear();

    this->numberOfVertices = numberOfVertices;
    for (int i = 0; i < numberOfVertices; i++)
    {
        this->treeStructure[i] = std::vector<int>();
    }
}

void GraphStructure::matchNodeByEdges(std::string &args)
{

    std::regex edgesRegex("([0-9]+,[0-9]+)");
    std::smatch edgeMatch;
    bool error = false;
    while (std::regex_search(args, edgeMatch, edgesRegex))
    {
        std::string aMatch = edgeMatch.str();
        // std::cout << aMatch << std::endl;
        std::stringstream aMatchStream(aMatch);
        std::string numberStr;
        std::vector<int> matchAsVector;

        while (!aMatchStream.eof())
        {
            std::getline(aMatchStream, numberStr, ',');
            int number = std::stoi(numberStr);

            if (this->treeStructure.find(number) == this->treeStructure.end())
            {
                std::cerr << "Error: Vertex " << number << " in {" << aMatch << "}"
                          << " is not found in graph" << std::endl;
                error = true;
                break;
            }
            matchAsVector.push_back(number);
        }

        if (error)
        {
            for (int i = 0; i < this->numberOfVertices; i++)
            {
                std::vector<int> newArr;
                this->treeStructure[i] = newArr;
            }
            break;
        }
        this->addToEdgeList(matchAsVector[0], matchAsVector[1]);
        this->treeStructure[matchAsVector[0]].push_back(matchAsVector[1]);
        this->treeStructure[matchAsVector[1]].push_back(matchAsVector[0]);

        args = edgeMatch.suffix();
    }
}

void GraphStructure::findPathBetweenDifferentNodes(int &fromNode, int &toNode, bool &found, std::vector<int> &resultPath)
{

    this->treeNodeLevel.clear();
    this->nodeQueue.push(fromNode);
    this->treeNodeLevel[fromNode] = 0;
    while (!this->nodeQueue.empty())
    {
        int currentSearch = this->nodeQueue.front();
        this->nodeQueue.pop();

        if (this->nodeSet.find(currentSearch) == this->nodeSet.end())
        {
            if (currentSearch == toNode)
            {
                found = true;
                // work node back to get path
                resultPath = this->findPath(currentSearch, fromNode, this->treeStructure, this->treeNodeLevel);
                break;
            }
            else
            {
                for (int i = 0; i < this->treeStructure[currentSearch].size(); i++)
                {
                    int nextSearch = this->treeStructure[currentSearch][i];
                    this->nodeQueue.push(nextSearch);
                    if (this->treeNodeLevel[nextSearch] == 0 && nextSearch != fromNode)
                    {
                        this->treeNodeLevel[nextSearch] = this->treeNodeLevel[currentSearch] + 1;
                    }
                }
                this->nodeSet.insert(currentSearch);
            }
        }
    }
}

void GraphStructure::findPathBetweenSameNodes(int &fromNode, int &toNode, bool &found, std::vector<int> &resultPath)
{
    std::set<std::string> checkedNeighbours;
    int minimumLength = this->treeStructure.size();
    for (int i = 0; i < this->treeStructure[fromNode].size(); i++)
    {
        for (int j = 0; j < this->treeStructure[fromNode].size(); j++)
        {
            if (i == j)
            {
                continue;
            }
            int neighbour1 = this->treeStructure[fromNode][i];
            int neighbour2 = this->treeStructure[fromNode][j];
            std::string pairNeighbours;

            if (neighbour1 < neighbour2)
            {
                pairNeighbours = std::to_string(neighbour1) + "," + std::to_string(neighbour2);
            }
            else
            {
                pairNeighbours = std::to_string(neighbour2) + "," + std::to_string(neighbour1);
            }

            if (checkedNeighbours.find(pairNeighbours) != checkedNeighbours.end())
            {
                continue;
            }
            bool tempFound = false;
            std::vector<int> tempResultPath;
            this->nodeSet.clear();
            this->nodeQueue = std::queue<int>();
            this->nodeSet.insert(fromNode);
            this->findPathBetweenDifferentNodes(neighbour1, neighbour2, tempFound, tempResultPath);

            if (tempFound)
            {
                if (tempResultPath.size() >= 2 && tempResultPath.size() < minimumLength)
                {
                    minimumLength = tempResultPath.size();
                    resultPath = tempResultPath;
                    found = tempFound;
                }
            }
        }
    }

    if (found)
    {
        resultPath.insert(resultPath.begin(), fromNode);
        resultPath.push_back(fromNode);
    }
}

std::vector<int> GraphStructure::findPath(int currentNode, int startingNode, std::map<int, std::vector<int>> treeStructure, std::map<int, int> treeNodeLevel)
{
    std::vector<int> reversePath;
    reversePath.push_back(currentNode);
    int currentLevel = treeNodeLevel[currentNode];
    std::vector<int> result;
    // std::cout << reversePath[0] << std::endl;

    do
    {
        int nextNode = treeStructure[currentNode][0];
        for (int i = 0; i < treeStructure[currentNode].size(); i++)
        {
            int neighbourNode = treeStructure[currentNode][i];
            if (treeNodeLevel.find(neighbourNode) != treeNodeLevel.end() && treeNodeLevel[neighbourNode] < currentLevel)
            {
                nextNode = neighbourNode;
                currentLevel = treeNodeLevel[neighbourNode];
            }
        };
        // std::cout << "Next node: " << nextNode << std::endl;
        reversePath.push_back(nextNode);
        currentNode = nextNode;
    } while (currentNode != startingNode);

    // std::cout << "length of reverse path is " << reversePath.size() << std::endl;
    for (int i = (reversePath.size() - 1); i >= 0; i--)
    {
        result.push_back(reversePath[i]);
    }

    return result;
}

bool getCommandsAndArgument(std::string line, char &command, std::string &args)
{
    std::istringstream input(line);

    std::ws(input);

    if (input.eof())
    {
        return false;
    }

    input >> command;

    std::ws(input);

    std::getline(input, args);

    return true;
}

void printMap(std::map<int, std::vector<int>> treeStructure)
{
    std::cout << "Tree structure : {" << std::endl;
    for (auto kv = treeStructure.begin(); kv != treeStructure.end(); ++kv)
    {
        std::cout << kv->first << ": "
                  << "[";
        for (auto v : kv->second)
        {
            std::cout << v << ", ";
        }
        std::cout << "]," << std::endl;
    }
    std::cout << "}" << std::endl;
}

void printTreeLevel(std::map<int, int> treeLevelNode)
{
    std::cout << "Node Level : {" << std::endl;
    for (auto kv = treeLevelNode.begin(); kv != treeLevelNode.end(); ++kv)
    {
        std::cout << kv->first << ": "
                  << kv->second << std::endl;
    }
    std::cout << "}" << std::endl;
}

class SatSolver
{
private:
    int max_vertices = 0;
    int starting_vertex = 1;
    std::vector<int> minimumVertexCover = std::vector<int>();
    std::unique_ptr<Minisat::Solver> solver;
    Minisat::vec<Minisat::Lit> clause;
    std::vector<std::vector<Minisat::Lit>> all_literals;
    std::vector<int> vertexEdges;
    void formula_1(int no_of_vertex, int k);
    void formula_1_v2(int no_of_vertex, int k);
    void formula_2(int no_of_vertex, int k);
    void formula_3(int no_of_vertex, int k);
    void formula_4(int k, std::vector<int> vertexEdges);
    void formula_4_v2(int k, std::vector<int> vertexEdges);

public:
    SatSolver(int new_max_vertices, std::vector<int> vertexEdges);
    std::vector<int> get_min_vertex_cover_cnf();
    std::vector<int> get_min_vertex_cover_3cnf();
};

SatSolver::SatSolver(int new_max_vertices, std::vector<int> new_vertexEdges)
{
    this->max_vertices = new_max_vertices;
    this->vertexEdges = new_vertexEdges;
    this->solver = std::unique_ptr<Minisat::Solver>(new Minisat::Solver());
}

void SatSolver::formula_1(int no_of_vertex, int currentMinimumNo)
{
    for (int i = 0; i < currentMinimumNo; i++)
    {
        for (int j = 0; j < no_of_vertex; j++)
        {
            this->clause.push(this->all_literals[j][i]);
        }
        this->solver->addClause(clause);
        clause.clear();
    }
}
void SatSolver::formula_1_v2(int no_of_vertex, int currentMinimumNo)
{

    for (int i = 0; i < currentMinimumNo; i++)
    {

        std::vector<Minisat::Lit> newVariables = std::vector<Minisat::Lit>(no_of_vertex - 1);
        for (int i = 0; i < no_of_vertex - 1; i++)
        {
            newVariables[i] = Minisat::mkLit(this->solver->newVar());
        }
        this->solver->addClause(this->all_literals[0][i], newVariables[0]);
        this->solver->addClause(this->all_literals[no_of_vertex - 1][i], ~newVariables[newVariables.size() - 1]);
        for (int j = 0; j < no_of_vertex - 2; j++)
        {
            this->clause.push(~newVariables[j]);
            this->clause.push(this->all_literals[j + 1][i]);
            this->clause.push(newVariables[j + 1]);
            this->solver->addClause(clause);
            clause.clear();
        }
    }
}

void SatSolver::formula_2(int no_of_vertex, int currentMinimumNo)
{
    for (int m = 0; m < no_of_vertex; m++)
    {
        for (int p = 0; p < currentMinimumNo - 1; p++)
        {
            for (int q = p + 1; q < currentMinimumNo; q++)
            {
                this->solver->addClause(~this->all_literals[m][p], ~this->all_literals[m][q]);
            }
        }
    }
}

void SatSolver::formula_3(int no_of_vertex, int currentMinimumNo)
{
    for (int m = 0; m < currentMinimumNo; m++)
    {
        for (int p = 0; p < no_of_vertex - 1; p++)
        {
            for (int q = p + 1; q < no_of_vertex; q++)
            {
                this->solver->addClause(~this->all_literals[p][m], ~this->all_literals[q][m]);
            }
        }
    }
}

void SatSolver::formula_4(int currentMinimumNo, std::vector<int> input_vertexEdges)
{
    for (int i = 0; i < vertexEdges.size(); i = i + 2)
    {
        for (int j = 0; j < currentMinimumNo; j++)
        {
            this->clause.push(this->all_literals[input_vertexEdges[i]][j]);
            this->clause.push(this->all_literals[input_vertexEdges[i + 1]][j]);
        }
        this->solver->addClause(clause);
        this->clause.clear();
    }
}

void SatSolver::formula_4_v2(int currentMinimumNo, std::vector<int> input_vertexEdges)
{

    for (int i = 0; i < vertexEdges.size(); i = i + 2)
    {
        std::vector<Minisat::Lit> newVariables = std::vector<Minisat::Lit>((2 * currentMinimumNo) - 1);
        for (int i = 0; i < newVariables.size(); i++)
        {
            newVariables[i] = Minisat::mkLit(this->solver->newVar());
        }
        int count = 0;
        for (int j = 0; j < currentMinimumNo; j++)
        {
            if (j == 0)
            {
                this->solver->addClause(newVariables[0], this->all_literals[input_vertexEdges[i]][j]);
                this->solver->addClause(~newVariables[newVariables.size() - 1], this->all_literals[input_vertexEdges[i + 1]][j]);
                continue;
            }
            this->clause.push(~newVariables[count]);
            this->clause.push(this->all_literals[input_vertexEdges[i]][j]);
            this->clause.push(newVariables[count + 1]);
            this->solver->addClause(clause);
            this->clause.clear();
            count += 1;
            this->clause.push(~newVariables[count]);
            this->clause.push(this->all_literals[input_vertexEdges[i + 1]][j]);
            this->clause.push(newVariables[count + 1]);
            this->solver->addClause(clause);
            this->clause.clear();
            count += 1;
        }
    }
}

std::vector<int> SatSolver::get_min_vertex_cover_cnf()
{
    if (this->vertexEdges.size() < 1)
    {
        return this->minimumVertexCover;
    }
    int max = this->max_vertices;

    for (int currentMinimumNo = 1; currentMinimumNo <= max; currentMinimumNo++)
    {
        // clear data of previous iterations
        this->minimumVertexCover.clear();
        this->clause.clear();
        this->solver.reset(new Minisat::Solver());
        this->all_literals = std::vector<std::vector<Minisat::Lit>>(max, std::vector<Minisat::Lit>(currentMinimumNo));

        // initialize literals
        for (int i = 0; i < max; i++)
        {
            for (int j = 0; j < currentMinimumNo; j++)
            {
                this->all_literals[i][j] = Minisat::mkLit(this->solver->newVar());
            }
        }

        // formula 1
        this->formula_1(max, currentMinimumNo);

        // formula 2
        this->formula_2(max, currentMinimumNo);

        // formula 3
        this->formula_3(max, currentMinimumNo);

        // formula 4
        this->formula_4(currentMinimumNo, this->vertexEdges);

        bool isSat = this->solver->solve();

        if (isSat == 1)
        {
            for (int i = 0; i < max; i++)
            {
                for (int j = 0; j < currentMinimumNo; j++)
                {
                    if (this->solver->modelValue(this->all_literals[i][j]) == (Minisat::lbool) true)
                    {
                        this->minimumVertexCover.push_back(i);
                    }
                }
            }
            std::sort(this->minimumVertexCover.begin(), this->minimumVertexCover.end());
            return this->minimumVertexCover;
        }
    }

    return this->minimumVertexCover;
}

std::vector<int> SatSolver::get_min_vertex_cover_3cnf()
{
    if (this->vertexEdges.size() < 1)
    {
        return this->minimumVertexCover;
    }
    int max = this->max_vertices;

    for (int currentMinimumNo = 1; currentMinimumNo <= max; currentMinimumNo++)
    {
        // clear data of previous iterations
        this->minimumVertexCover.clear();
        this->clause.clear();
        this->solver.reset(new Minisat::Solver());
        this->all_literals = std::vector<std::vector<Minisat::Lit>>(max, std::vector<Minisat::Lit>(currentMinimumNo));

        // initialize literals
        for (int i = 0; i < max; i++)
        {
            for (int j = 0; j < currentMinimumNo; j++)
            {
                this->all_literals[i][j] = Minisat::mkLit(this->solver->newVar());
            }
        }

        // formula 1
        this->formula_1_v2(max, currentMinimumNo);

        // formula 2
        this->formula_2(max, currentMinimumNo);

        // formula 3
        this->formula_3(max, currentMinimumNo);

        // formula 4
        this->formula_4_v2(currentMinimumNo, this->vertexEdges);

        bool isSat = this->solver->solve();

        if (isSat == 1)
        {
            for (int i = 0; i < max; i++)
            {
                for (int j = 0; j < currentMinimumNo; j++)
                {
                    if (this->solver->modelValue(this->all_literals[i][j]) == (Minisat::lbool) true)
                    {
                        this->minimumVertexCover.push_back(i);
                    }
                }
            }
            std::sort(this->minimumVertexCover.begin(), this->minimumVertexCover.end());
            return this->minimumVertexCover;
        }
    }

    return this->minimumVertexCover;
}

int main(int argc, char **argv)
{
    // Test code. Replaced with your code
    // new
    GraphStructure graphStructure;
    graphStructure.numberOfVertices = 0;

    do
    {

        // new
        graphStructure.treeNodeLevel = std::map<int, int>();
        graphStructure.nodeQueue = std::queue<int>();
        graphStructure.nodeSet = std::set<int>();
        graphStructure.edgeList.clear();

        std::string line;
        char command;
        std::string args;
        std::getline(std::cin, line);

        if (!getCommandsAndArgument(line, command, args))
        {
            continue;
        }

        if (command == 'V')
        {

            // new
            graphStructure.generateGraphMap(std::stoi(args));
            // printMap(treeStructure);
        }
        else if (command == 'E')
        {
            // new
            graphStructure.matchNodeByEdges(args);
            // printMap(treeStructure);
            SatSolver newSatSolver(graphStructure.numberOfVertices, graphStructure.edgeList);

            // get mininum vertex cover
            auto result = newSatSolver.get_min_vertex_cover_3cnf();

            for (int i = 0; i < result.size(); i++)
            {
                std::cout << result[i] << " ";
            }

            std::cout << std::endl;
        }
        else if (command == 's')
        {
            // new
            std::string resultPathStr = "";
            std::vector<int> resultPath;
            std::istringstream fromAndToStream(args);
            int fromNode;
            int toNode;
            fromAndToStream >> fromNode;
            std::ws(fromAndToStream);
            fromAndToStream >> toNode;
            bool found = false;

            if (graphStructure.treeStructure.find(fromNode) == graphStructure.treeStructure.end())
            {
                std::cerr << "Error: Vertex " << fromNode << " does not exist in graph" << std::endl;
            }
            else if (graphStructure.treeStructure.find(toNode) == graphStructure.treeStructure.end())
            {
                std::cerr << "Error: Vertex " << toNode << " does not exist in graph" << std::endl;
            }
            else
            {
                if (toNode != fromNode)
                {

                    graphStructure.findPathBetweenDifferentNodes(fromNode, toNode, found, resultPath);
                }
                else
                {
                    graphStructure.findPathBetweenSameNodes(fromNode, toNode, found, resultPath);
                }
                if (!found)
                {
                    std::cerr << "Error: could not find path betweeen " << fromNode << " and " << toNode << std::endl;
                }
                else
                {
                    for (int i = 0; i < resultPath.size() - 1; i++)
                    {
                        resultPathStr += std::to_string(resultPath[i]) + "-";
                    }
                    resultPathStr += std::to_string(resultPath[resultPath.size() - 1]);
                    std::cout << resultPathStr << std::endl;
                }
            }
        }
        else
        {
            std::cerr << "Error: not a valid command" << std::endl;
        }
    } while (!std::cin.eof());

    return 0;
}
