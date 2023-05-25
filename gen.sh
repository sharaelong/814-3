clang++ -o generator -std=c++17 -O2 generator.cpp 

for i in {1..50}; do ./generator > "data/data$i.in"; done
