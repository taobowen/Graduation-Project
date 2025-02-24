# Design of Ant Algorithm for Solving Wait-Free Pipeline Scheduling Problem
## Wait-Free Pipeline Scheduling Problem
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm1.jpg)
1) A machine can only process one workpiece at a time.
2) Each workpiece is processed at most once on each machine.
3) Each machine processes workpieces in the same order.
4) The processing time of each workpiece on a given machine is fixed.
For all such and other scheduling problems based on this, the research goal is to find a specific sequence of workpieces that minimizes the total processing time of the sequence.
## Ant Algorithm
The ant algorithm is essentially a probabilistic algorithm for solving combinatorial optimization problems. Its bionic basis comes from the path search behavior of ants during foraging. The algorithm has the characteristics of parallel computing, heuristic search, and positive information feedback. The ant algorithm is widely used in TCP problems to solve the optimal route of the traveling salesman. It takes advantage of its probabilistic algorithm characteristics of finding the optimal path. The search for the optimal workpiece sequence in the no-wait assembly line scheduling problem can also be quantified as a path search process. Starting from the first workpiece, the ant selects the next workpiece in turn based on information such as pheromones and heuristic functions. According to the calculation formula for the difference in completion time between two given workpieces in the sequence, we can derive a "path length" for any given sequence. Taking the total processing time as the objective function value, we can find the "shortest path".
## Orthogonal test
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm2.jpg)
## Algorithm Overview
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm3.jpg)
