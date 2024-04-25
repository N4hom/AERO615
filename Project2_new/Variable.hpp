#include <iostream>
#include <vector>
#include <string>

template<typename T>
class Variable {
private:
    std::string name_;
    std::vector<T> data_;

public:
    Variable(const std::string& name) : name_(name) {}

    // Custom resize method
    void resize(size_t size) {
        data_.resize(size);
    }

    // Custom resize method for 2D vectors
    void resize(size_t rows, size_t cols) {
        data_.resize(rows, std::vector<T>(cols));
    }

    // Accessors
    T& operator[](size_t index) {
        return data_[index];
    }

    const T& operator[](size_t index) const {
        return data_[index];
    }

    // Getter for name
    std::string name() const {
        return name_;
    }
};
