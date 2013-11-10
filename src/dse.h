#ifndef _dse_h_
#define _dse_h_

#include <stdint.h>
#include <map>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <queue>
#include <Eigen/SparseCore>
#include <parallel/algorithm>
#include "DisjointSet.h"

/* Entry from the matrix M */
typedef struct {
    uint64_t i;
    uint64_t j;
    double value;
} C_tuple;

/* How we compare two tuples */
bool C_tuple_compare(const C_tuple *c1, const C_tuple *c2)
{
    bool retVal = false;

    if (c1->value < c2->value)
    {
        retVal = true;
    }
    else if (c1->value == c2->value)
    {
        if (c1->i < c2->i)
        {
            retVal = true;
        }
        else if (c1->i == c2->i)
        {
            if (c1->j < c2->j)
            {
                retVal = true;
            }
        }
    }

    return retVal;
}

/* Node in the tree for building the dendrogram */
typedef struct node {
    node* parent;
    node* leftChild;
    node* rightChild;
    uint64_t vertex;
    uint64_t numEdges;
    uint64_t numNodes;
    double density;
} node;

/* Read the graph file into a sparse matrix */
Eigen::SparseMatrix<double> readGraph(std::string filename)
{
    std::ifstream graph;
    graph.open(filename);

    if (!graph)
    {
        std::cerr << "Error reading graph" << std::endl;
        exit(1);
    }

    uint64_t n, m;
    graph >> n >> m;

    Eigen::SparseMatrix<double> mat(n, n) ;
    typedef Eigen::Triplet<double> Ty;
    std::vector<Ty> tripletList;
    tripletList.reserve(2*m);

    uint64_t v1, v2;
    for (uint64_t i = 0; i < m; i++)
    {
        graph >> v1 >> v2;
        tripletList.push_back(Ty(v1,v2,1));
        tripletList.push_back(Ty(v2,v1,1));
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

/* Find all of the entries in the matrix M and store them in the C tuple */
void getC(Eigen::SparseMatrix<double> &g, std::vector<C_tuple*> &C)
{
    /* Normalize every column of X */
    auto X(g);
    for (int32_t i = 0; i < X.cols(); i++)
    {
        X.col(i) /= X.col(i).norm();
    }

    /* Compute M */
    Eigen::SparseMatrix<double> Z = X.transpose();
    Eigen::SparseMatrix<double> M = (Z * X).triangularView<Eigen::Upper>();

    /* Set diagonal to zeros */
    for (int32_t i = 0; i < X.rows(); i++)
        M.coeffRef(i, i) = 0;

    /* Get all C_tuples from M */
    for (int32_t i = 0; i < M.outerSize(); ++i)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it)
        {
            double val = it.value();
            if (val != 0)
            {
                C_tuple *c = new C_tuple;
                c->i = it.row();
                c->j = it.col();
                c->value = val;
                C.push_back(c);
            }
        }
    }
}

/* Print a node to stdout */
void printNode(node* n)
{
    std::cout << "vertex: " << n->vertex << "\t#nodes: " << n->numNodes << "\t#edges: " << n->numEdges << "\tdensity: " << n->density << std::endl;
}

/* Do a post order traversal of the tree rooted at node root, and perform
 * the visit function on each node */
void postorder(void (*visit)(node*), node* root)
{
    if (root->leftChild)  postorder(visit, root->leftChild);
    if (root->rightChild) postorder(visit, root->rightChild);
    visit(root);
}

/* Perform a level order traversal and visit each node */
void levelorder(void (*visit)(node*), node* root)
{
    std::queue<node*> q;
    q.push(root);

    node* n;
    while(!q.empty())
    {
        n = q.front();
        q.pop();
        if (n->leftChild)  q.push(n->leftChild);
        if (n->rightChild) q.push(n->rightChild);
        visit(n);
    }
}

/* Create the tree for the dendrogram representation */
node* createT(Eigen::SparseMatrix<double> &g, std::vector<C_tuple*> C)
{
    /* Nodes that represent the actual vertices in the graph */
    node *vertices = new node[g.rows()];

    for (int32_t i = 0; i < g.rows(); i++)
    {
        vertices[i].parent = NULL;
        vertices[i].leftChild = NULL;
        vertices[i].rightChild = NULL;
        vertices[i].vertex = i;
        vertices[i].numNodes = 1;
        vertices[i].numEdges = 0;
        vertices[i].density = 0.0;
    }

    /* Nodes that will be created as part of the construction of T */
    node **upperTree = new node*[g.rows()];
    for (int32_t i = 0; i < g.rows(); i++)
        upperTree[i] = &(vertices[i]);

    /* DisjointSet object with n elements */
    DisjointSet d(g.rows());

    for (auto it = C.crbegin(); it != C.crend(); ++it)
    {
        C_tuple *entry = *it;
        uint64_t setI = d.findSet(entry->i);
        uint64_t setJ = d.findSet(entry->j);

        if (setI == setJ) continue;
        else d.unionSets(setI, setJ);

        node *v1 = upperTree[setI];
        node *v2 = upperTree[setJ];

        node *n = new node;
        n->parent = NULL;
        n->vertex = -1;     // This evaluates to roughly 4 billion
        n->numNodes = 0;
        n->numEdges = 0;
        n->density = 0.0;

        if (v1->numNodes >= v2->numNodes)
        {
            n->leftChild = v1;
            n->rightChild = v2;
        }
        else
        {
            n->leftChild = v2;
            n->rightChild = v1;
        }

        v1->parent = n;
        v2->parent = n;

        uint32_t sn = d.findSet(entry->i);
        upperTree[sn] = n;
    }

    //if (!root) std::cerr << "root is null" << std::endl;
    node* root = NULL;//vertices[0].parent;
    int32_t i = 0;
    while (i < g.rows() && !root) root = vertices[i++].parent;
    if (i >= g.rows())
    {
        std::cerr << "Shit. Something went wrong and we fail." << std::endl;
        exit(1);
    }
    while (root->parent) root = root->parent;

    /* This root is the root node for the tree T, which is the dendrogram hierarchy. */
    return root;
}

/* Find the lowest common ancestor for nodes n1 and n2 */
node* LCA(node* n1, node* n2)
{
    if (!n1)
    {
        std::cerr << "returning null" << std::endl;
        return NULL;
    }
    if (!n2)
    {
        std::cerr << "returning null" << std::endl;
        return NULL;
    }

    std::vector<node*> p1;
    std::vector<node*> p2;

    /* Build the root to n1 path */
    node* temp = n1;
    while(temp->parent)
    {
        p1.push_back(temp->parent);
        temp = temp->parent;
    }

    /* Build the root to n2 path */
    temp = n2;
    while(temp->parent)
    {
        p2.push_back(temp->parent);
        temp = temp->parent;
    }

    // handle the case when one of them is the root
    if (p1.size() == 0)
        return n1;
    else if (p2.size() == 0)
        return n2;


    /* Find last spot where they have the same value */
    int64_t i = p1.size()-1, j = p2.size()-1;
    node *lca = NULL;
    while (p1[i] == p2[j])
    {
        if (i == 0)
        {
            lca = p1[i];
            break;
        }
        else if (j == 0)
        {
            lca = p2[j];
            break;
        }
        else if (p1[i-1] != p2[j-1])
        {
            lca = p1[i];
            break;
        }

        i--; j--;
    }

    return lca;
}

static void addToMap(std::map<uint64_t, node*> &m, node* n)
{
    if (n->leftChild) addToMap(m, n->leftChild);
    if (n->rightChild) addToMap(m, n->rightChild);
    if (n->vertex != (uint64_t)-1) m[n->vertex] = n;
}

static void countVerticesAndEdgesWrapUp(node* root)
{
    if (root->leftChild && root->rightChild)
    {
        countVerticesAndEdgesWrapUp(root->leftChild);
        countVerticesAndEdgesWrapUp(root->rightChild);
        root->numEdges = root->leftChild->numEdges + root->rightChild->numEdges + root->numEdges;
        root->numNodes = root->leftChild->numNodes + root->rightChild->numNodes;
    }
    else
    {
        root->numNodes = 1;
    }
}

/* Count how many edges and vertices each subtree has */
void countVerticesAndEdges(Eigen::SparseMatrix<double> &g, node* root)
{
    std::map<uint64_t, node*> v2n;
    addToMap(v2n, root);

    postorder([] (node *n) {n->numEdges = 0;}, root);


    Eigen::SparseMatrix<double> gU = g.triangularView<Eigen::Upper>();
    //#pragma omp parallel for schedule(dynamic)
    for (int32_t i = 0; i < gU.outerSize(); ++i)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(gU, i); it; ++it)
        {
            uint64_t i = it.row();
            uint64_t j = it.col();

            node* lca = LCA(v2n[i], v2n[j]);
            if (lca) lca->numEdges++;
        }
    }

    countVerticesAndEdgesWrapUp(root);
}

/* Compute the density of a node */
void computeDensity(node* root)
{
    if (root->leftChild) computeDensity(root->leftChild);
    if (root->rightChild) computeDensity(root->rightChild);
    if (root->numNodes == 1) root->density = 0;
    else root->density = (((double)root->numEdges) / (((double)(root->numNodes * (root->numNodes-1))) / 2.0));
}

/* Extract all the subgraphs that are above the given threshold recursively */
void extractSubgraphs(node* root, double t)
{
    if (root->density > t)
    {
        postorder([] (node *n) {if (!n->leftChild && !n->rightChild) std::cout << n->vertex << " ";}, root);
        std::cout << std::endl;
    }
    else if (root->leftChild && root->rightChild)
    {
        extractSubgraphs(root->leftChild, t);
        extractSubgraphs(root->rightChild, t);
    }
}

#endif /* _dse_h_ */
