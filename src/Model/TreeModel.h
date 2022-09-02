/*
 * TreeModel.h
 *
 *  Created on: August 17, 2022
 *      Author: hideakiv
 */

#ifndef TREEMODEL_H_
#define TREEMODEL_H_

#include "Node.h"
#include <map>

template <class T> 
class TreeModel {
public:
    std::map<int,Node<T>*> nodes;
public:
    TreeModel();
    ~TreeModel() { }
    void print_graph();
    void add_node(Node<T>& n);
    void add_child(T obj, Node<T>& p);
};
template <class T> void TreeModel<T>::print_graph() {
    typename std::map<int,Node<T>*>::iterator itr;
    for (itr = nodes.begin(); itr!=nodes.end(); ++itr) {
        itr->second->print_node();
    }
}
template <class T> void TreeModel<T>::add_node(Node<T>& n) {
    int node_id = n.id();
    nodes[node_id] = &n;
}
template <class T> void TreeModel<T>::add_child(T obj, Node<T>& p) {
    int node_id = nodes.size();
    Node<T> c(node_id,p,obj);
    nodes[node_id] = &c;
}

#endif //TREEMODEL_H_