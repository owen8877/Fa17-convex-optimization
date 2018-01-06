#include <iostream>
#include <string>

using namespace std::string_literals;

struct A {
    A()
    : history("created")
    {
    }

    A(A&& r)
    : history("move-constructed,"s + r.history)
    {
        r.history = "zombie: was "s + r.history;
    }
    A(const A& r)
    : history("copied from: " + r.history)
    {
    }
    ~A() {
        history = "destroyed,"s + history;
        std::cout << history << std::endl;
    }
    A& operator=(A&& r) {
        history = "move-assigned from " + r.history + " (was "s + history + ")"s;
        r.history = "zombie: was "s + r.history;
        return *this;
    }
    A& operator=(const A&r ) {
        history = "copied from " + r.history;
        return *this;
    }
    std::string history;
};

A foo() {
    auto p = std::make_pair(A{}, 2);
    // ... do something
    return p.first;
}



auto main() -> int
{
    auto a = foo();
    return 0;
}