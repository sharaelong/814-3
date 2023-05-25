# Progress

1. 2-approximation and square grid  
Score: 36515513.01  
Time: 12ms, relatively fast.  
It is well-known so chatgpt do it well. Also I assigned first 19 cities to extra salesmen, and divide whole grid 11x11 square and assign each square to each saleman.  

2. Improve square grid to use 10x14 grid  
Score: 32342713.38  
Time: 12ms, same as before  
It uses every salesmen in tour.  

3. 2-opt heuristic  
Score: 26611309.87  
Time: 552ms  
I applied new technique to find more improved TSP, no change to clustering. ChatGPT gives a code about it. Also I made 32ms time limit for each salesman instance so it will not exceed time limit. However real time is 552ms which means every subproblem was solved under 32ms.  

4. Clustering with K-means  
Score: 34219073.69  
Time: 1120ms  
I used k-means for clustering (which also given by ChatGPT), but it shows worse results than naive 10x14 grid. Max iteration was 100, and clustering takes time about 550ms purely. I don't have any intuition about how k-means groups points in random distribution...  

5. 2-opt with random initial tour
Score: 26519866.54  
Time: 760ms  
I called 2-opt function multiple time with random initial tour. But solution quality is not improved at much.  

6. 3-opt heuristic  
Score: 25835278.68  
Time: 1384ms  
I used 3-opt heuristic instead of multiple 2-opt heuristic. It showed 2.6% performance gain, which is similar as other resource's experiment results. Check every 3-opt moves and select best move. However, I can see execution time is increased sharply.  

# Utilities
- `score.sh`: compile main source code and score calculator, print score. `tmp-file` is temporary file for saving result.
- `gen.sh`: re-generate test cases.

# Useful resources
- [Traveling salesman problem (wikipedia)](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
- [TSP Basics](http://tsp-basics.blogspot.com)
- [The Traveling Salesman Problem: A Computational Study, David L. Applegate, Robert E. Bixby, Vasek Chvátal & William J. Cook](https://www.math.uwaterloo.ca/tsp/book/index.html)
- [Aarts E., Lenstra J.K. Local Search in combinatorial optimization](https://www.amazon.com/Local-Search-Combinatorial-Optimization-Emile/dp/0691115222)
- [Combining 2-opt, 3-opt and 4-opt with k-swap-kick perturbations for the traveling salesman problem](https://isd.ktu.lt/it2011/material/Proceedings/1_AI_5.pdf)
