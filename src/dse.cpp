#include <iostream>
#include "dse.h"
#include <math.h>
#include <parallel/algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Core>

#define DEBUG 0

int main(int argc, char *argv[])
{
    std::cerr << "reading graph...";
    Eigen::SparseMatrix<double> g = readGraph(std::string(argv[1]));
    std::cerr << "done." << std::endl;

#if DEBUG
    std::cout << "Adjacency Matrix" << std::endl;
    std::cout << g << std::endl;
#endif

    std::vector<C_tuple*> C;
    std::cerr << "getting C vector...";
    getC(g, C);

    std::cerr << "sorting...";
    __gnu_parallel::sort(C.begin(), C.end(), C_tuple_compare);
    //std::sort(C.begin(), C.end(), C_tuple_compare);
    std::cerr << "done." << std::endl;

#if DEBUG
    for (uint64_t i = 0; i < C.size(); i++)
    {
        std::cout << C[i]->i << ", " << C[i]->j << ", " << C[i]->value<< std::endl;
    }
#endif

    std::cerr << "creating T...";
    node* root = createT(g, C);
    std::cerr << "done." << std::endl;

#if DEBUG
    //postorder(printNode, root);
    levelorder(printNode, root);

    node* left = root->leftChild;
    while (left->leftChild) left = left->leftChild;

    node* right = root->rightChild;
    while (right->rightChild) right = right->rightChild;

    node* lca = LCA(left, right);
    std::cout << "lca for " << left->vertex << ", " << right->vertex << ": " << lca->vertex << std::endl;

    left = right->parent->leftChild;
    lca = LCA(left, right);
    std::cout << "lca for " << left->vertex << ", " << right->vertex << ": " << lca->vertex << std::endl;
#endif

    std::cerr << "counting vertices and edges...";
    countVerticesAndEdges(g, root);
    std::cerr << "done." << std::endl;
    std::cerr << "computing density...";
    computeDensity(root);
    std::cerr << "done." << std::endl;

#if DEBUG
    std::cout << "\nPrinting after the countVerticesAndEdges\n" << std::endl;
    postorder(printNode, root);
#endif

    std::cerr << "extracting subgraphs...";
    extractSubgraphs(root, 0.75);
    std::cerr << "done." << std::endl;
}
