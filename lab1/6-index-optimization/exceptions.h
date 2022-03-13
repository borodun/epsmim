#ifndef LAB1_EXCEPTIONS_H
#define LAB1_EXCEPTIONS_H

#include <exception>
#include <string>

class WaveSimException : public std::exception {
private:
    std::string errorString;
public:
    explicit WaveSimException(std::string errStr) {
        errorString = std::move(errStr);
    }

    const char *what() const noexcept override {
        return errorString.c_str();
    }
};

class badArgsException : public WaveSimException {
public:
    explicit badArgsException(std::string errStr) : WaveSimException(std::move(errStr)) {}
};

class fileException : public WaveSimException {
public:
    explicit fileException(std::string errStr) : WaveSimException(std::move(errStr)) {}
};

#endif
