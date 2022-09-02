/*
 * Node.h
 *
 *  Created on: August 17, 2022
 *      Author: hideakiv
 */

#ifndef NODE_H_
#define NODE_H_


#include <map>
#include <iostream>


template <class T> 
class Node {
private:
    const int id_num;
public:
    Node<T>* parent;
    std::map<int,Node<T>*> children;
    T info;

public:
    Node(int id, Node<T>& p, T& obj):
        id_num(id),
        parent(&p),
        info(obj)
    { }
    ~Node(){ }
    bool operator == (const Node &n) const { return id_num == n.id(); }
    bool operator != (const Node &n) const { return !operator==(n); }
    int id() const { return id_num; }
    void add_child(Node<T>& n);

    void print_node();
};

template <class T> void Node<T>::add_child(Node<T>& n) {
    int id = n.id();
    children[id] = &n;
}

template <class T> void Node<T>::print_node() {
    std::cout << "node ID: " << id_num <<"\n";
    std::cout << "\tparent ID: " << parent->id() <<"\n";
    std::cout << "\tchild  ID: ";
    typename std::map<int,Node<T>*>::iterator itr;
    for (itr = children.begin(); itr!=children.end(); itr++) {
        std::cout << itr->second->id() << ", ";
    }
    std::cout << "\n";
}

#endif //NODE_H_