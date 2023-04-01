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

std::vector<int> APPROX_VC_1(std::map<int, std::vector<int>> edge_dic)
{
    std::vector<int> cover;
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
    return cover;
}

std::vector<int> REFINED_APPROX_VC_1(std::map<int, std::vector<int>> edge_dic_2)
{
    std::vector<int> new_cover;
    std::vector<int> cover = APPROX_VC_1(edge_dic_2);
    std::vector<int> cover_2 = cover;

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
    return new_cover;
}

std::vector<int> APPROX_VC_2(std::map<int, std::vector<int>> edge_dic)
{
    std::vector<int> cover;
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
    return cover;
}

std::vector<int> REFINED_APPROX_VC_2(std::map<int, std::vector<int>> edge_dic_2)
{
    std::vector<int> new_cover;
    std::vector<int> cover = APPROX_VC_2(edge_dic_2);
    std::vector<int> cover_2 = cover;

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
    return new_cover;
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
                        std::cout << "CNF_vertex_cover: ";
                        std::vector<int> result = vertex_cover(edge_values, Limit); // cnf implementation
                        sort(result.begin(), result.end());

                        if (result.size() > 0)
                        {
                            for (int i = 0; i < result.size(); i++)
                            {
                                std::cout << result[i] << " ";
                            }
                            std::cout << std::endl;
                        }

                        std::vector<int> cover = APPROX_VC_1(edge_dic);

                        std::cout << "APPROX-VC-1: ";
                        for (const auto &v : cover)
                        {
                            std::cout << v << " ";
                        }

                        std::cout << "\n";

                        std::vector<int> new_cover = REFINED_APPROX_VC_1(edge_dic);
                        std::cout << "REFINED-APPROX-VC-1: ";
                        for (const auto &v : new_cover)
                        {
                            std::cout << v << " ";
                        }

                        std::cout << "\n";
                        std::cout << "\n";

                        std::vector<int> cover_1 = APPROX_VC_2(edge_dic);
                        std::cout << "APPROX-VC-2: ";
                        for (const auto &v : cover_1)
                        {
                            std::cout << v << " ";
                        }

                        std::vector<int> new_cover_1 = REFINED_APPROX_VC_2(edge_dic);
                        std::cout << "\n";
                        std::cout << "REFINED-APPROX-VC-2: ";
                        for (const auto &v : new_cover_1)
                        {
                            std::cout << v << " ";
                        }
                        std::cout << "\n";
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