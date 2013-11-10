#ifndef _DISJOINT_SET_H_
#define _DISJOINT_SET_H_

#include <cstddef>
#include <stdint.h>
#include <vector>

class DisjointSet
{
    private:
        /* Private node class */
        struct Node
        {
            uint32_t index;     // The index this node represents. For us this
                                // will be the vertex ID.
            uint32_t rank;      // The rank of this node
            Node* parent;       // The parent of this node
        };

        /* Vector to store all of the elements in */
        std::vector<Node*> elements;

        /* Total number of sets */
        uint32_t nSets;

    public:
        DisjointSet();
        DisjointSet(uint32_t nElements);    // Create a DisjointSet with n elements initially
        ~DisjointSet();

        /* Returns the total number of elements. */
        uint32_t numElements() const;

        /* Returns the number of sets. */
        uint32_t numSets() const;

        /* Find the set that this element is associated with. */
        uint32_t findSet(uint32_t) const;

        /* Add a new element, it will be in its own set. */
        void addElement();

        /* Add N new elements, each in their own set */
        void addElementBatch(uint32_t);

        /* Union two sets into one. */
        void unionSets(uint32_t, uint32_t);
};

#endif /* _DISJOINT_SET_H_ */
