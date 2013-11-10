#include "DisjointSet.h"

DisjointSet::DisjointSet()
{
    nSets = 0;
}

DisjointSet::DisjointSet(uint32_t n)
{
    nSets = 0;
    addElementBatch(n);
}

DisjointSet::~DisjointSet()
{
    for (uint32_t i = 0; i < elements.size(); i++)
    {
        delete elements[i];
    }
}

uint32_t DisjointSet::numElements() const
{
    return elements.size();
}

uint32_t DisjointSet::numSets() const
{
    return nSets;
}

uint32_t DisjointSet::findSet(uint32_t e) const
{
    if (e >= elements.size()) throw "Element out of range in findSet.";

    // Find the root of the set e is in
    Node* curr = elements[e];
    while (curr->parent) curr = curr->parent;

    // Now that we know the root for the set, update all the other nodes to
    // point directly to it so that in the future a query can be answered
    // quicker.
    Node* root = curr;
    Node* next = NULL;
    curr = elements[e];
    while (curr != root)
    {
        next = curr->parent;
        curr->parent = root;
        curr = next;
    }

    return root->index;
}

void DisjointSet::addElement()
{
    Node *n = new Node;
    n->index = elements.size();
    n->rank = 0;
    n->parent = NULL;
    elements.push_back(n);
}

void DisjointSet::addElementBatch(uint32_t n)
{
    // Current number of elements
    uint32_t nElements = elements.size();

    // insert a batch of place holder elements for now
    elements.insert(elements.end(), n, (Node*) NULL);

    // fill out the placeholders
    for (uint32_t i = nElements; i < elements.size(); i++)
    {
        Node *n = new Node;
        n->index = i;
        n->rank = 0;
        n->parent = NULL;
        elements[i] = n;
    }

    // We now have n more sets since each new element is in their own set
    nSets += n;
}

void DisjointSet::unionSets(uint32_t s1, uint32_t s2)
{
    if (s1 >= elements.size() || s2 >= elements.size())
        throw "Set out of range in unionSets.";

    Node* set1 = elements[s1];
    Node* set2 = elements[s2];

    if (set1->rank > set2->rank)
    {
        set2->parent = set1;
    }
    else if (set1->rank < set2->rank)
    {
        set1->parent = set2;
    }
    else
    {
        set1->parent = set2;
        set1->rank++;
    }

    nSets--;
}
