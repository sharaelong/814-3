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

7. Dynamic clustering  
Score: 25241052.9  
Time: 4352ms  
This clustering algorithm tries moving points from a randomly selected cluster to a randomly selected adjacent cluster. Cluster is initialized with 10x14 grid, so neighbors of certain cluster refer to those connected on four sides. The selection probability is exponential value to current cluster's tour length. When point moving increases max tour length between two cluster, rollback that changes.  
Some lessons I got from this:
- BOJ is much slower than my local machine. It's easy to overestimate algorithm performance.
- Just moving single point is really weak.
- Exponential probability of choice is good heuristic, indeed it seems necessary.
- 3-opt is not good to use as a blackbox for its huge variance.
- Heuristic depends on area of cluster convex hull is not good, because sharp shape has small area but increases tour length significantly.

8. Some fancy dynamic clustering  
Score: Catastrophe!  
Time: Also too long!  
I used many approaches in this step. But everything failed. Inscribed circle area of delaunay triangulation, MST length of group, connectivity constraint in every move are not enough to avoid points floated over space wide. Inherently triangulation become weird at the edge of the region. To meet this issue I tried to define 'adjacency' between two points at maximum distance of d, but it was not good approach either.  

9. Cut edges of maximum length cluster  
Score: 25454746.62  
Time: 3956ms  
When initialize cluster with 10x14 grid, its tour length shows kind of like normal distribution with sigma about 30000. So we can think of cut the edges of two biggest cluster and distribute edges to adjacent grid will be good at lower maximum tour length.  

10. Shift rectangular grid
Score: 24768752.05  
Time: 3856ms  
Simple 10x14 grid gives approximate normal distribution on tour length in each cluster. So removing outlier gives big performance improvement. From this observation, this algorithm traverse every row in grid. In each row find biggest and smallest cluster by tour length, and move fragments length of `DX` between that clusters. Other clusters between that clusters are 'shifted', preserving their area. Also same process will operated for every column, shifting fragments length of `DY`.  

11. Shift rectangular grid (Advanced)
Score: 24604028.2  
Time: 3924ms  
This is a extension of trial number 10. First in y-axis (14 cuts), calculate every tour length of initial grid and accumulate tour length which in same y-axis grid. Then regard as 14 metal bars with different mass. Cut it into 14 new bars with same mass. After that do the same thing in x-axis, but this time each y-line is independent. We have to do cutting only inside the same y-line.
p.s. o3 and ofast don't bring us fast execution time. Maybe it's becuse reversing array is the most bottleneck of our code.

FYI) I embedded a [LKH solver](http://akira.ruc.dk/~keld/research/LKH/). For now `Gain23` function is weird in my embedding, so I just turned off. I have to test some options to optimize performance without losing ability to get optimum solution.

12. Embed LKH Solver
Score: 23241597.44  
Time: 3000ms  
I finally successed embedding LKH solver, and it gives me about 6% of performance improvements. It seemed to give optimum solution in many instances, which leads my shifting algorithm can believe tour value. I used one level tree, maxcand=4, patchingA=1, patchingC=0 option with 1 run, 1 trial. Also it is worth to note that solver is even faster than my 3-opt code.

# Remaining Ideas
- Grouping some points before any calculation. For example we can assume that two points within a distance of 1000 are expected to belong in same group at last. Then treat this grouped points as single point, 7th approach will become better because single point movement is actually multiple points movement.

# Utilities
- `score.sh`: compile main source code and score calculator, print score. `tmp-file` is temporary file for saving result.
- `gen.sh`: re-generate test cases.
- `visualizer.py`: visualize cluster result.

# Useful resources
- [Traveling salesman problem (wikipedia)](https://en.wikipedia.org/wiki/Travelling_salesman_problem)
- [TSP Basics](http://tsp-basics.blogspot.com)
- [The Traveling Salesman Problem: A Computational Study, David L. Applegate, Robert E. Bixby, Vasek Chvátal & William J. Cook](https://www.math.uwaterloo.ca/tsp/book/index.html)
- [Aarts E., Lenstra J.K. Local Search in combinatorial optimization](https://www.amazon.com/Local-Search-Combinatorial-Optimization-Emile/dp/0691115222)
- [Combining 2-opt, 3-opt and 4-opt with k-swap-kick perturbations for the traveling salesman problem](https://isd.ktu.lt/it2011/material/Proceedings/1_AI_5.pdf)
