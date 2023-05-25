* Progress

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
