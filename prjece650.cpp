#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <regex>
#include <set>
#include <queue>
#include <tuple>
#include <memory>
#include "minisat/core/SolverTypes.h"
#include "minisat/core/Solver.h"
#include <pthread.h>
#include <signal.h>
#include <unistd.h>

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

struct CNF_Vals
{
    int noVertices;
    std::vector<int> edgList;
    bool foundSolution;
    std::vector<int> result;
};

struct vc_Vals
{
    std::map<int, std::vector<int>> edgList;
    std::vector<int> result;
};

void *CNF_THREAD(void *args)
{
    struct CNF_Vals *data = (struct CNF_Vals *)args;
    // Install signal handler for SIGALRM
    struct sigaction sa;
    sa.sa_handler = [](int sig)
    { pthread_exit(NULL); };
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGALRM, &sa, NULL);
    SatSolver newSatSolverCNF(data->noVertices, data->edgList);
    // Set the alarm for 30 seconds
    alarm(30);
    data->result = newSatSolverCNF.get_min_vertex_cover_cnf();
    data->foundSolution = true;
    // Disable the alarm
    alarm(0);
    return NULL;
}

void *CNF3_THREAD(void *args)
{
    struct CNF_Vals *data = (struct CNF_Vals *)args;
    // Install signal handler for SIGALRM
    struct sigaction sa;
    sa.sa_handler = [](int sig)
    { pthread_exit(NULL); };
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGALRM, &sa, NULL);
    SatSolver newSatSolverCNF(data->noVertices, data->edgList);
    // Set the alarm for 30 seconds
    alarm(30);
    data->result = newSatSolverCNF.get_min_vertex_cover_3cnf();
    data->foundSolution = true;
    // Disable the alarm
    alarm(0);
    return NULL;
}
// std::vector<int> APPROX_VC_1(std::map<int, std::vector<int>> edge_dic)
void *APPROX_VC_1(void *args)
{

    std::vector<int> cover;
    vc_Vals *data = (vc_Vals *)args;
    std::map<int, std::vector<int>> edge_dic = (*data).edgList;
    while (edge_dic.size() != 0)
    {
        int max_key = 0;
        int max_len = 0;

        for (const auto &p : edge_dic)
        {
            int key = p.first;
            int val = p.second.size();
            if (val > max_len)
            {
                max_key = key;
                max_len = val;
            }
        }

        cover.push_back(max_key);

        for (const auto &value : edge_dic[max_key])
        {
            auto it = std::find(edge_dic[value].begin(), edge_dic[value].end(), max_key);
            if (it != edge_dic[value].end())
            {
                edge_dic[value].erase(it);
            }
            if (edge_dic[value].size() == 0)
            {
                edge_dic.erase(value);
            }
        }
        edge_dic.erase(max_key);
    }
    sort(cover.begin(), cover.end());
    data->result = cover;
    return NULL;
}

// std::vector<int> REFINED_APPROX_VC_1(std::map<int, std::vector<int>> edge_dic_2)
void *REFINED_APPROX_VC_1(void *args)
{
    std::vector<int> new_cover;
    APPROX_VC_1(args);
    vc_Vals *data = (vc_Vals *)args;
    std::map<int, std::vector<int>> edge_dic_2 = (*data).edgList;

    std::vector<int> cover = std::vector<int>((*data).result);
    ;
    std::vector<int> cover_2 = std::vector<int>((*data).result);
    ;

    if (cover.size() == 1)
    {
        new_cover.push_back(cover[0]);
    }
    else
    {
        for (auto x : cover_2)
        {
            bool should_add_x = false;
            for (auto y : edge_dic_2[x])
            {
                if (std::find(cover.begin(), cover.end(), y) == cover.end())
                {
                    should_add_x = true;
                    break;
                }
            }
            if (should_add_x)
            {
                new_cover.push_back(x);
            }
            else
            {
                cover.erase(std::remove(cover.begin(), cover.end(), x), cover.end());
            }
        }
    }
    sort(new_cover.begin(), new_cover.end());
    data->result = new_cover;
    return NULL;
}

// std::vector<int> APPROX_VC_2(std::map<int, std::vector<int>> edge_dic)
void *APPROX_VC_2(void *args)
{
    std::vector<int> cover;
    vc_Vals *data = (vc_Vals *)args;
    std::map<int, std::vector<int>> edge_dic = (*data).edgList;
    while (!edge_dic.empty())
    {
        int first_key = edge_dic.begin()->first;
        int first_value = edge_dic.begin()->second[0];
        cover.push_back(first_key);
        cover.push_back(first_value);
        for (auto &value : edge_dic[first_key])
        {
            auto &other_values = edge_dic[value];
            other_values.erase(std::remove(other_values.begin(), other_values.end(), first_key), other_values.end());
            if (other_values.empty())
            {
                edge_dic.erase(value);
            }
        }
        edge_dic.erase(first_key);
        for (auto &value : edge_dic[first_value])
        {
            auto &other_values = edge_dic[value];
            other_values.erase(std::remove(other_values.begin(), other_values.end(), first_value), other_values.end());
            if (other_values.empty())
            {
                edge_dic.erase(value);
            }
        }
        edge_dic.erase(first_value);
    }
    sort(cover.begin(), cover.end());
    data->result = cover;
    return NULL;
}

// std::vector<int> REFINED_APPROX_VC_2(std::map<int, std::vector<int>> edge_dic_2)
void *REFINED_APPROX_VC_2(void *args)
{
    std::vector<int> new_cover;
    APPROX_VC_2(args);
    vc_Vals *data = (vc_Vals *)args;
    std::map<int, std::vector<int>> edge_dic_2 = (*data).edgList;

    std::vector<int> cover = std::vector<int>((*data).result);
    std::vector<int> cover_2 = std::vector<int>((*data).result);
    if (cover.size() == 1)
    {
        new_cover.push_back(cover[0]);
    }
    else
    {
        for (auto x : cover_2)
        {
            bool should_add_x = false;
            for (auto y : edge_dic_2[x])
            {
                if (std::find(cover.begin(), cover.end(), y) == cover.end())
                {
                    should_add_x = true;
                    break;
                }
            }
            if (should_add_x)
            {
                new_cover.push_back(x);
            }
            else
            {
                cover.erase(std::remove(cover.begin(), cover.end(), x), cover.end());
            }
        }
    }
    sort(new_cover.begin(), new_cover.end());
    data->result = new_cover;
    return NULL;
}

void mini_sat_clauses_1(std::vector<int> edge_values, int Limit, std::vector<std::vector<Minisat::Lit>> &vector_of_vector, Minisat::Solver *solver)
{
    for (int i = 0; i < vector_of_vector[i].size(); i++)
    {
        Minisat::vec<Minisat::Lit> sub_vec;
        for (int j = 0; j < vector_of_vector.size(); j++)
        {
            sub_vec.push(vector_of_vector[j][i]);
        }
        solver->addClause(sub_vec);
    }
}

void mini_sat_clauses_2(std::vector<int> edge_values, int Limit, std::vector<std::vector<Minisat::Lit>> &vector_of_vector, Minisat::Solver *solver)
{
    for (int i = 0; i < vector_of_vector.size(); i++)
    {
        for (int j = 0; j < vector_of_vector[i].size() - 1; j++)
        {
            for (int k = j + 1; k < vector_of_vector[k].size(); k++)
            {
                solver->addClause(~vector_of_vector[i][j], ~vector_of_vector[i][k]);
            }
        }
    }
}

void mini_sat_clauses_3(std::vector<int> edge_values, int Limit, std::vector<std::vector<Minisat::Lit>> &vector_of_vector, Minisat::Solver *solver)
{
    for (int i = 0; i < vector_of_vector[i].size(); i++)
    {
        for (int j = 0; j < Limit - 1; j++)
        {
            for (int k = j + 1; k < vector_of_vector.size(); k++)
            {
                solver->addClause(~vector_of_vector[j][i], ~vector_of_vector[k][i]);
            }
        }
    }
}

void mini_sat_clauses_4(std::vector<int> edge_values, int Limit, std::vector<std::vector<Minisat::Lit>> &vector_of_vector, Minisat::Solver *solver)
{
    for (int i = 0; i < edge_values.size(); i += 2)
    {
        Minisat::vec<Minisat::Lit> sub_vec;
        for (int j = 0; j < vector_of_vector[j].size(); j++)
        {
            sub_vec.push(vector_of_vector[edge_values[i]][j]);
            sub_vec.push(vector_of_vector[edge_values[i + 1]][j]);
        }
        solver->addClause(sub_vec);
    }
}

void mini_sat_clauses_1_4(std::vector<int> edge_values, int Limit, std::vector<std::vector<Minisat::Lit>> vector_of_vector, Minisat::Solver *solver)
{
    mini_sat_clauses_1(edge_values, Limit, vector_of_vector, solver);
    mini_sat_clauses_2(edge_values, Limit, vector_of_vector, solver);
    mini_sat_clauses_3(edge_values, Limit, vector_of_vector, solver);
    mini_sat_clauses_4(edge_values, Limit, vector_of_vector, solver);
}

std::vector<int> vertex_cover(std::vector<int> edge_values, int Limit)
{
    std::vector<int> vertex_cover_1;

    for (int cover_length = 1; cover_length <= Limit; cover_length++)
    {
        std::unique_ptr<Minisat::Solver> solver(new Minisat::Solver());
        std::vector<std::vector<Minisat::Lit>> vector_of_vector(Limit);

        for (int i = 0; i < Limit; i++)
        {
            for (int j = 0; j < cover_length; j++)
            {
                vector_of_vector[i].push_back(Minisat::mkLit(solver->newVar()));
            }
        }

        mini_sat_clauses_1_4(edge_values, Limit, vector_of_vector, solver.get());

        if (solver->solve())
        {
            for (int i = 0; i < vector_of_vector.size(); i++)
            {
                for (int j = 0; j < vector_of_vector[j].size(); j++)
                {
                    if (solver->modelValue(vector_of_vector[i][j]) == (Minisat::lbool) true)
                    {
                        vertex_cover_1.push_back(i);
                    }
                }
            }
            return vertex_cover_1;
        }
        solver.reset(new Minisat::Solver());
    }
}

std::vector<int> short_part(const std::map<int, std::vector<int>> edge_dic, int start, int end)
{
    if (start == end)
    {
        std::vector<std::vector<int>> path_list;
        int range_1 = (edge_dic.at(end)).size() - 1;
        for (int i = 0; i < range_1; i++)
        {
            int right = edge_dic.at(end)[i];
            int range_2 = (edge_dic.at(end)).size();
            for (int j = i + 1; j < range_2; j++)
            {
                int left = edge_dic.at(end)[j];
                std::queue<std::vector<int>> queue;
                std::vector<int> path;
                std::map<int, bool> checked_val;
                checked_val[end] = true;
                checked_val[right] = true;
                queue.push({right});
                while (!queue.empty())
                {
                    path = queue.front();
                    queue.pop();
                    int end_node_in_queue = path.back();
                    if (end_node_in_queue == left)
                    {
                        path_list.push_back(path);
                    }
                    for (const auto element : edge_dic.at(end_node_in_queue))
                    {
                        if (checked_val[element])
                        {
                            continue;
                        } // Check viited elements
                        checked_val[element] = true;
                        std::vector<int> new_path = path;
                        new_path.push_back(element);
                        queue.push(new_path);
                    }
                }
            }
        }
        int path_list_size = path_list.size();
        if (path_list_size > 0)
        {
            int smallest_vector_size = path_list[0].size();
            for (const auto vectors_within : path_list)
            {
                int current_size = vectors_within.size();
                if (current_size < smallest_vector_size)
                {
                    smallest_vector_size = vectors_within.size();
                }
            }
            for (const auto vectors_within : path_list)
            {
                int current_size_2 = vectors_within.size();
                if (current_size_2 == smallest_vector_size)
                {
                    return vectors_within;
                }
            }
        }
        return {};
    }
    else
    {
        std::queue<std::vector<int>> queue;
        std::vector<int> path;
        std::map<int, bool> checked_val;
        checked_val[start] = true;
        queue.push({start});
        while (!queue.empty())
        {
            path = queue.front();
            queue.pop();
            int end_node_in_queue = path.back();
            if (end_node_in_queue == end)
            {
                return path;
            }
            for (const auto element : edge_dic.at(end_node_in_queue))
            {
                if (checked_val[element])
                {
                    continue;
                }
                checked_val[element] = true;
                std::vector<int> new_path = path;
                new_path.push_back(element);
                queue.push(new_path);
            }
        }
        return {};
    }
}

int main()
{
    int Limit = 0;
    std::set<int> set_keys;
    std::map<int, std::vector<int>> edge_dic;
    std::vector<int> edge_values;

    while (!std::cin.eof())
    {

        std::string line;
        std::getline(std::cin, line);
        std::istringstream input(line);
        std::ws(input);
        std::string command;
        input >> command;

        if (command == "V")
        {
            std::ws(input);
            std::string arguments;
            input >> arguments;
            Limit = std::stoi(arguments);
            edge_dic = {};
            set_keys = {};
            edge_values = {};
            // std::cout<<Limit<<std::endl;
        }
        else if (command == "E")
        {
            if (Limit > 0)
            {
                std::ws(input);
                std::string arguments;
                input >> arguments;
                arguments = arguments.substr(1, arguments.size() - 2);
                std::vector<std::tuple<int, int>> key_val;
                std::regex parse_line("<(\\d+),(\\d+)>");
                auto start_line = std::sregex_iterator(arguments.begin(), arguments.end(), parse_line);
                auto end_line = std::sregex_iterator();
                for (std::sregex_iterator i = start_line; i != end_line; ++i)
                {
                    std::smatch match = *i;
                    int x = std::stoi(match[1].str());
                    int y = std::stoi(match[2].str());
                    key_val.push_back({x, y});
                }
                int count = 0;
                for (const auto val : key_val)
                {
                    int x, y;
                    std::tie(x, y) = val;
                    if (x >= Limit || y >= Limit)
                    {
                        // std::cout <<x<<" "<<y<<" = "<< Limit<<"\n";
                        count++;
                    }
                };

                if (count > 0)
                {
                    std::cout << "Error: An Edge value is greater than the provided limit, Provide a valid Edge Case" << std::endl;
                }
                else
                {
                    for (const auto val : key_val)
                    {
                        int x, y;
                        std::tie(x, y) = val;
                        edge_dic[x].push_back(y);
                        edge_dic[y].push_back(x);
                        set_keys.insert(x);
                        set_keys.insert(y);
                    };

                    for (const auto val : key_val)
                    {
                        int x, y;
                        std::tie(x, y) = val;
                        edge_values.push_back(x);
                        edge_values.push_back(y);
                    };
                    /*
                    for (int i=0;i<edge_values.size();i++){
                        std::cout<<edge_values[i]<<std::endl;
                    }
                    */
                    if (set_keys.size() > 0)
                    {
                        // thread id
                        pthread_t CNF_SAT, CNF_3SAT, VC1, VC2, VC1_REF, VC2_REF;

                        // create thread

                        // CNF_SAT CREATE THREAD
                        // SatSolver newSatSolverCNF(Limit, edge_values);
                        CNF_Vals cnf_args = {.noVertices = Limit, .edgList = edge_values, .foundSolution = false, .result = std::vector<int>()};
                        pthread_create(&CNF_SAT, NULL, &CNF_THREAD, &cnf_args);

                        // CNF_3SAT CREATE THREAD
                        CNF_Vals cnf3_args = {.noVertices = Limit, .edgList = edge_values, .foundSolution = false};
                        pthread_create(&CNF_3SAT, NULL, &CNF3_THREAD, &cnf3_args);

                        // VC1 CREATE THREAD
                        vc_Vals vc_1_arg = {.edgList = edge_dic, .result = std::vector<int>()};
                        // std::map<int, std::vector<int>> edge_dic_APPROX_VC_1(edge_dic);
                        pthread_create(&VC1, NULL, &APPROX_VC_1, &vc_1_arg);
                        // VC2 CREATE THREAD
                        vc_Vals vc_2_arg = {.edgList = edge_dic, .result = std::vector<int>()};
                        // std::map<int, std::vector<int>> edge_dic_APPROX_VC_2(edge_dic);
                        pthread_create(&VC2, NULL, &APPROX_VC_2, &vc_2_arg);
                        // VC1_REF CREATE THREAD
                        vc_Vals vc_1_ref_arg = {.edgList = edge_dic, .result = std::vector<int>()};
                        // std::map<int, std::vector<int>> edge_dic_REFINED_APPROX_VC_1(edge_dic);
                        pthread_create(&VC1_REF, NULL, &REFINED_APPROX_VC_1, &vc_1_ref_arg);
                        // VC2_REF CREATE THREAD
                        vc_Vals vc_2_ref_arg = {.edgList = edge_dic, .result = std::vector<int>()};
                        // std::map<int, std::vector<int>> edge_dic_REFINED_APPROX_VC_2(edge_dic);
                        pthread_create(&VC2_REF, NULL, &REFINED_APPROX_VC_2, &vc_2_ref_arg);

                        // join thread

                        pthread_join(CNF_SAT, NULL);
                        pthread_join(CNF_3SAT, NULL);
                        pthread_join(VC1, NULL);
                        pthread_join(VC2, NULL);
                        pthread_join(VC1_REF, NULL);
                        pthread_join(VC2_REF, NULL);

                        std::vector<int> cover = vc_1_arg.result;
                        std::vector<int> cover_1 = vc_2_arg.result;
                        std::vector<int> new_cover = vc_1_ref_arg.result;
                        std::vector<int> new_cover_1 = vc_2_ref_arg.result;
                        // print the output of each thread

                        // CNF PRINT
                        std::cout << "CNF-SAT-VC: ";
                        // std::vector<int> result_cnf = newSatSolverCNF.get_min_vertex_cover_cnf(); // cnf implementation
                        std::vector<int> result_cnf = cnf_args.result;
                        if (cnf_args.foundSolution == false)
                        {
                            std::cout << "timeout";
                        }
                        else
                        {
                            if (result_cnf.size() > 0)
                            {
                                std::cout << result_cnf[0];
                                for (int i = 1; i < result_cnf.size(); i++)
                                {
                                    std::cout << "," << result_cnf[i];
                                }
                            }
                        }
                        std::cout << std::endl;

                        // CNF-3 PRINT
                        std::cout << "CNF-3-SAT-VC: ";
                        // std::vector<int> result_cnf = newSatSolverCNF.get_min_vertex_cover_cnf(); // cnf implementation
                        std::vector<int> result_3cnf = cnf3_args.result;
                        if (cnf_args.foundSolution == false)
                        {
                            std::cout << "timeout";
                        }
                        else
                        {
                            if (result_3cnf.size() > 0)
                            {
                                std::cout << result_3cnf[0];
                                for (int i = 1; i < result_3cnf.size(); i++)
                                {
                                    std::cout << "," << result_3cnf[i];
                                }
                            }
                        }
                        std::cout << std::endl;
                        // VC-1 PRINT
                        // std::vector<int> cover = APPROX_VC_1(edge_dic);
                        std::cout << "APPROX-VC-1: ";
                        if (cover.size() > 0)
                        {
                            std::cout << cover[0];
                            for (int i = 1; i < cover.size(); i++)
                            {
                                std::cout << "," << cover[i];
                            }
                        }
                        std::cout << std::endl;
                        // VC-2 PRINT
                        // std::vector<int> cover_1 = APPROX_VC_2(edge_dic);
                        std::cout << "APPROX-VC-2: ";
                        if (cover_1.size() > 0)
                        {
                            std::cout << cover_1[0];
                            for (int i = 1; i < cover_1.size(); i++)
                            {
                                std::cout << "," << cover_1[i];
                            }
                        }
                        std::cout << std::endl;

                        // VC-1-REF PRINT
                        // std::vector<int> new_cover = REFINED_APPROX_VC_1(edge_dic);
                        std::cout << "REFINED-APPROX-VC-1: ";
                        if (new_cover.size() > 0)
                        {
                            std::cout << new_cover[0];
                            for (int i = 1; i < new_cover.size(); i++)
                            {
                                std::cout << "," << new_cover[i];
                            }
                        }
                        std::cout << std::endl;

                        // VC-2-REF PRINT
                        // std::vector<int> new_cover_1 = REFINED_APPROX_VC_2(edge_dic);
                        std::cout << "REFINED-APPROX-VC-2: ";
                        if (new_cover_1.size() > 0)
                        {
                            std::cout << new_cover_1[0];
                            for (int i = 1; i < new_cover_1.size(); i++)
                            {
                                std::cout << "," << new_cover_1[i];
                            }
                        }
                        std::cout << std::endl;
                    }
                    else
                    {
                        std::cout << std::endl;
                    }
                }
            }
            else
            {
                std::cout << "Error: Please provide V limit" << std::endl;
            }
        }
        else if (command == "s")
        {
            std::ws(input);
            std::string arguments_1;
            input >> arguments_1;
            int start = std::stoi(arguments_1);
            std::ws(input);
            std::string arguments_2;
            input >> arguments_2;
            int end = std::stoi(arguments_2);
            if (set_keys.size() > 0)
            {
                if (set_keys.count(start) > 0 && set_keys.count(end) > 0)
                {
                    std::vector<int> path = short_part(edge_dic, start, end);

                    if (path.empty())
                    {
                        std::cout << "Error: There is no path from the start node to the end node." << std::endl;
                    }

                    else if (start == end)
                    {
                        std::cout << start << "-";
                        for (int val : path)
                        {
                            std::cout << val;
                            if (val != path[path.size() - 1])
                            {
                                std::cout << "-";
                            }
                        }
                        std::cout << "-" << end;
                        std::cout << "\n";
                    }
                    else
                    {
                        for (int val : path)
                        {
                            std::cout << val;
                            if (val != path[path.size() - 1])
                            {
                                std::cout << "-";
                            }
                        }
                        std::cout << "\n";
                    }
                }
                else
                {
                    std::cout << "Error: One or both Vertex Value provided doesn't exist" << std::endl;
                }
            }
            else
            {
                std::cout << "Error: There are no edge values provided" << std::endl;
            }
        }
    }
}