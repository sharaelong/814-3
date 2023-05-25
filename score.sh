clang++ -o main -std=c++17 -O2 main.cpp
clang++ -o score -std=c++17 -O2 score.cpp

rm tmp-file
touch tmp-file
for i in {1..50}; do cat "data/data$i.in" | ./main | ./score "$i" >> tmp-file || break; done
paste -sd+ tmp-file | sed 's/+$//g' | bc

