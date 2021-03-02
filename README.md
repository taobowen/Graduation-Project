# 求解无等待流水线调度问题的蚂蚁算法设计
## 无等待流水线调度问题
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm1.jpg)
1)	一台机器最多只能同时加工一个工件。
2)	每个工件在每台机器上最多只会加工一次。
3)	每台机器加工工件的顺序相同。
4)	每个工件在给定机器上的加工时间固定。
对于所有此类以及以此为基础的调度问题，研究目标都是找到一个特定的工件序列，使得序列的总加工时间最短。
## 蚂蚁算法
蚂蚁算法本质上是一种用于解决组合优化问题的概率型算法。它的仿生学基础源自于蚂蚁觅食过程中的路径搜索行为。该算法具有并行计算、启发式搜索和信息正反馈的特点。蚂蚁算法被大量应用于TCP问题去求解旅行商的最优路线正是利用了其寻找优化路径的概率型算法的特点，而无等待流水线调度问题中最优工件序列的搜索也可以量化为一种路径搜索的过程。从第一个工件开始，蚂蚁根据信息素、启发函数等信息依次选择下一个工件，根据序列中给定两工件之间的完工时间之差的计算公式，我们可以得出任意给定序列的一个“路径长度”，以总加工时间为目标函数值，就可以求出“最短路径”。
## 正交测试
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm2.jpg)
## 算法总览
![avatar](https://github-picture.oss-cn-beijing.aliyuncs.com/ant_algorithm3.jpg)
