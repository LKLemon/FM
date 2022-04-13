#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <string>

int main() {
   std::vector<bool> a = {1,1,0,1};
   std::vector<bool> b = {0,1,1,0};
   std::vector<bool> c(4,1);
   std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<bool>());
    c.push_back(a.back());
    c.back() = 0;

    std::cout << a.back() << "\n";

    return 0;
}
