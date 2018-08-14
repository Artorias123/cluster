cluster
模式识别实验，聚类程序

聚类方法包括：
K-Means、Fuzzy C-Means、Gaussian Mixed Model、Fast Search & Find Density Peaks

使用方法：
将待分类数据放入clusterdata中
在set.txt中设置参数，运行run
在result中即可查看聚类结果

Fast Search & Find Density Peaks算法开启后将先利用改算法得到初始聚类中心的坐标，再使用聚类方法进行聚类。
