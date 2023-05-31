# clang++ -o main -std=c++17 -O2 main.cpp
# clang++ -o score -std=c++17 -O2 score.cpp

# rm tmp-file
# touch tmp-file
# for i in {1..50}; do (cat "data/data$i.in" | ./main | ./score "$i" >> tmp-file) || break; done
# paste -sd+ tmp-file | sed 's/+$//g' | bc

clang++ -o main -std=c++17 -O2 main.cpp 
clang++ -o score -std=c++17 -O2 score.cpp

rm tmp-file
touch tmp-file

total_files=50
for ((i=1; i<=total_files; i++)); do
  printf '\rProcessing data%s.in ...' $i
  (cat "data/data$i.in" | ./main | ./score "$i" >> tmp-file) || break
done

if [[ $? -eq 0 ]]; then
  printf '\rCalculating score ...\n'
  paste -sd+ tmp-file | sed 's/+$//g' | bc
fi
