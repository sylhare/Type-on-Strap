---
published: true
tag: cpp
---
## Example

C++ have rich support manipulators for display.

[manipulators_cpp](https://www.cplusplus.com/reference/ios/noshowpos/)


```
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;


int main() {
double A = 100;
double B = 2006.0123;
double C = 23314159;

//LINE 1

cout << hex << left << showbase << nouppercase; // formatting
cout << (long long) A << endl; // actual printed part

// LINE 2
cout << dec << right << setw(15) << setfill('_') << showpos << fixed << setprecision(2); // formatting
cout << B << endl; // actual printed part

// LINE 3
cout << scientific << uppercase << noshowpos << setprecision(9); // formatting
cout << C << endl; // actual printed part

return 0;
}
```

output:  
```
0x64             
_______+2006.01  
2.331415927E+03
```
refer 1 question and answer on hackerank. [refer this exercise](https://www.hackerrank.com/challenges/prettyprint/)