#ifndef CORE_DATA_STRUCTURES_H
#define CORE_DATA_STRUCTURES_H

#include "memory.hpp"
#include <cstdlib>
#include <cstdio>
#include <type_traits>

#include <initializer_list>

namespace core {
namespace containers {

template<typename T>
bool order_default(const T& a, const T& b) {
    return a < b;
}

template<typename T>
class collection {
    T* data;
    unsigned size;
    unsigned reserved;

    void partition(bool(*compare)(const T& a, const T& b), int high, int low);
    void quicksort(bool(*compare)(const T& a, const T& b), int high, int low);
    void swap(int i, int j);

    unsigned index_of(const T& e);
    unsigned index_of(const T&& e);
public:
    static const unsigned out_of_range = -1;

    collection();
    collection(std::initializer_list<T> l);
    ~collection();

    unsigned count();

    bool have(const T& e);
    bool have(const T&& e);

    void order(bool(*compare)(const T& a, const T& b) = order_default<T>);
    
    void insert(const T e);

    void remove(const T&& e);
    void remove(T& e);

    void remove_index(unsigned index);
};

template<typename T>
collection<T>::collection() {
    this->size = 0;
    this->reserved = 10;
    this->data = (T*)memory::alloc(sizeof(T)*this->reserved);
}

template<typename T>
collection<T>::collection(std::initializer_list<T> l) {
    this->size = 0;
    this->reserved = 10;
    this->data = (T*)memory::alloc(sizeof(T)*this->reserved);
    
    for(auto i = l.begin(); i!= l.end(); i++)
        this->insert(*i);
}

template<typename T>
collection<T>::~collection() {
    if(std::is_pointer<T>::value)
        for(int i=0;i<this->size; i++)
            this->remove_index(0);
    
    free(this->data);
}

template<typename T>
unsigned collection<T>::count() {
    return this->size;
}

template<typename T>
bool collection<T>::have(const T& e) {
    return this->index_of(e);
}

template<typename T>
bool collection<T>::have(const T&& e) {
    return this->index_of(e);
}

template<typename T>
void collection<T>::swap(int i, int j) {
    T tmp = this->data[i];
    this->data[i] = this->data[j];
    this->data[j] = tmp;
}

template<typename T>
void collection<T>::partition(bool(*compare)(const T& a, const T& b), int high, int low) {
    int pivot = this->data[high];
    int i = low - 1;
    
    for(int j=low; j<high; j++)
        if(compare(this->data[j], this->data[pivot]))
            this->swap(++i, j);

    this->swap(i+1, high);
}

template<typename T>
void collection<T>::quicksort(bool(*compare)(const T& a, const T& b), int high, int low) {
    if(low < high) {
        int pi = this->partition(compare, low, high);
        quicksort(low, pi - 1);
        quicksort(pi, high);
    }
}

template<typename T>
void collection<T>::order(bool(*compare)(const T& a, const T& b)) {
    quicksort(compare, 0, this->size - 1);
}

template<typename T>
void collection<T>::insert(const T e) {
    if(this->reserved < this->size + 1) {
        this->reserved += size / 2;
        T* old = this->data;
        this->data = (T*)memory::alloc(sizeof(T)*this->reserved);
        memory::copy(this->data, old, this->size*sizeof(T));
        memory::free(old);
    }

    this->data[size++] = e;
}


template<typename T>
unsigned collection<T>::index_of(const T& e) {
    unsigned it = 0;

    while(it < this->size) {
        if(e == data[it])
            return it;

        it++;
    }

    return collection<T>::out_of_range;
}

template<typename T>
unsigned collection<T>::index_of(const T&& e) {
    unsigned it = 0;
    while(it < this->size) {
        if(std::is_pointer<T>::value && *e == *data[it]) {
            memory::free(e);
            return it;
        } else if(e == data[it]) {
            return it;
        }
        it++;
    }

    return collection<T>::out_of_range;
}

template<typename T>
void collection<T>::remove(const T&& e) {
    unsigned index = index_of(e);

    if(index == collection<T>::out_of_range)
        return;

    if(std::is_pointer<T>::value)
        memory::free(e);

    remove(index);
}

template<typename T>
void collection<T>::remove_index(unsigned index) {

    if(std::is_pointer<T>::value)
        memory::free(this->data[index]);

    this->data[index] = this->data[this->size];

    this->size--;

    unsigned diff = this->size - this->reserved;
    if(diff > 50) {
        this->reserved -= 50;
        T* old = this->data;
        
        this->data = memory::alloc(sizeof(T)*reserved - this->reserved);
        
        memory::copy(this->data, old, this->size);
        memory::free(old);
    }

}

template<typename T>
void collection<T>::remove(T& e) {
    unsigned index = index_of(e);

    if(index == collection<T>::out_of_range)
        return;

    remove_index(index);
}

} // namespace memory
} // namespace core

#endif