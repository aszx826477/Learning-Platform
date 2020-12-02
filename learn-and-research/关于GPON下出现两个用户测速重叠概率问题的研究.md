# 关于GPON下出现两个用户测试重叠概率问题的研究

## 1、问题描述

某GPON端口下有n=25（或n=28）个用户，每个用户每天测试的概率为p=0.5，每次测试持续时间为t=3分钟或者t=0.5分钟。集中测试时长在4个小时或者6个小时, 即T=240分钟或T=360分钟。求每天出现两个用户进行测试重叠的概率。

两个用户进行测试重叠定义：至少存在两个用户在测试时间上有重叠。

## 2、问题抽象

我们可以将问题分解和抽象为以下两个随机事件的概率问题：

* **问题一**：总共有$n$个用户，每个用户被抽取的概率为$p$，它们被抽取是相互独立的，求有$k(0 \leq k \leq n)$个用户被抽取到的概率

* **问题二**：在区间$[0, T]$上随机取$k$个点，求存在两点距离小于等于$t$的概率

### 2.1 求解问题一

问题一显然是一个二项分布（*Binomial Distribution*），我们将事件$A$记为“有用户被抽到”，随机变量$X$为事件$A$发生的次数，事件$A^k$表示事件A发生的次数为$k$，那么

$k$个用户被抽取到的概率服从二项分布，记为$X \sim B(n,p)$

$$P(A^k)=C_n^k \cdot p^k \cdot (1-p)^{n-k}$$

对于问题二，我们将事件$B$记为在区间$[0,T]$上随机取点存在两点距离小于等于$t$，事件$B^k$表示在在区间$[0,T]$上随机取$k$个点，存在两点距离小于等于$t$，则我们可以将问题二记为$P(B^k)$

显然，当$k=0,1$时，$P(B^k)=0$

那么每天出现两个用户在测试时间上有重叠的概率可以表示为

$$\sum_{k=2}^nP(A^k) \cdot P(B^k)$$

其中$P(A^k)$的表达式我们已经求出，那么问题的关键就在于求解$P(B^k)$

### 2.2 求解问题二

记$f_k(T,t) \triangleq P(\overline{B^k})$，所求$P(B^k)=1-f_k(T,t)$，即考虑在区间$[0,T]$上$k$个用户测速不存在重叠的概率是一个关于$T$和$t$的函数

已知

$$f_0(T,t)=f_1(T,t)=1$$

当$k=2$时，我们采用线性规划来求$P(B^2)$，如图所示，横轴和纵轴分别为事件一名测试用户$x$，另一名测试用户$y$开始的时间。当两者开始时间满足$|x-y| \leq t$的时候，两个测试时间会重叠。重叠部分即为图中阴影部分

![线性规划](pic/两个用户测速重叠的研究/线性规划.png)

根据几何概型，可求得

$$P(B^2)=1-\frac{(T-t)^2}{T^2}$$

$$f_2(T,t)=1-P(B^2)=1-[1-(\frac{T-t}{T})^2]=(1-\frac{t}{T})^2 \tag{2.1}$$

且易知

$$f_k(T,t)=0，对于k \geq [\frac{T}{t}]+1$$

对于一般的$k:2 < k \leq [\frac{T}{t}]$，记$x_i \sim U(0,T)(i=1,2,\dots,k), T_k=\underset{1 \leq i \leq k}{max} \ x_i$

$$f_k(T,t)=P(\overline{B^k})=\int_0^T P(\overline{B^k}|T_k=x) \cdot \rho_{T_k}(x)dx$$

其中$\rho_{T_k}(x)$为$T_k$的p.d.f.（概率密度函数）且

$$\rho_{T_k}(x)=\frac{kx^{k-1}}{T^k},0 \leq x \leq T$$

以及

$$P(\overline{B^k}|T_k=x)=0，对于0 \leq x \leq (k-1)t$$

且对$(k-1)t < x \leq T$，

$$
P(\overline{B^k}|T_k=x)
=P(k-1个用户在区间[0,x-t]上不存在重叠)
$$

$$
=f_{k-1}(x-t,t)
$$

可推出

$$
f_k(T,t)=\int_{(k-1)t}^T P(\overline{B^k}|T_k=x) \cdot \rho_{T_k}(x) dx \\
$$

$$
=\int_{(k-1)t}^T \frac{kx^{k-1}}{T^k} \cdot f_{k-1}(x-t,t)dx \tag{2.2}
$$

所以，我们可以得到关于$f_k(T,t)$的递推表达式

当$k=3$时，根据式(2.1)，可得

$$
f_3(T,t)
=\int_{2t}^T \frac{3x^2}{T^3} \cdot f_2(x-t,t)dx
$$

$$
=\int_{2t}^T \frac{3x^2}{T^3} \cdot (1-\frac{t}{x-t})^2 dx
$$

$$
=\frac{(2t-T)^3 \cdot (2t+T)}{(t-T) \cdot T^3} \tag{2.3}
$$


则

$$
P(B^3) = 1-f_3(T,t) = 1-\frac{(2t-T)^3 \cdot (2t+T)}{(t-T) \cdot T^3}
$$

当$k=4$时，根据式(2.3)，可得

$$
f_4(T,t)
=\int_{3t}^T \frac{4x^3}{T^4} \cdot f_3(x-t,t)dx
$$

$$
=\int_{3t}^T \frac{4x^3}{T^4} \cdot \frac{(2t-(x-t))^3  (2t+(x-t))}{(t-(x-t))  (x-t)^3} dx
$$


以此类推，根据$f_{k-1}(T,t)$我们可求得$f_k(T,t)$。而求出了$f_k(T,t)$即可求出$P(B^k)$

### 2.3 根据递推表达式求通项公式

求通项公式还不知道怎么做，只能利用递推公式求数值解，可以用蒙特卡洛模拟去求近似的数值解与理论解对比

## 3、蒙特卡洛模拟实验

对于每一个$k$使用python进行1万次的模拟实验

```
import random

total_num = 10000 # 模拟次数为10,000次
n = 25 # 一个GPON端口下接入的用户数量
p = 0.5 # 每个用户每天测试的概率
t = 3 # 每次测试的持续时间（分钟）
T = 240 # 测试总时长（分钟）


def simulate_one_time_k_users(k, t, T):
    user_list = list()

    for i in range(0, k):
        user_list.append(random.uniform(0, T))

    for i in range(len(user_list)):
        for j in range(i+1, len(user_list)):
            if abs(user_list[i] - user_list[j]) <= t:
                return 1 # 说明测试时间出现重叠

    return 0 # 说明没有出现测试时间重叠

# 某一天有k个用户进行测试
for k in range(2, n+1): # 某一天测试的用户数量
    overlap_num = 0 # 重叠的次数
    for i in range(total_num): # 进行多次的模拟实验，次数为total_num
        overlap_num += simulate_one_time_k_users(k, t, T)
    print("k = %s, T = %s, t = %s, overlap_prob = %.1f%%" % (k, T, t, overlap_num / total_num * 100))
```

模拟结果如下，`overlap_prob`是出现两个用户重叠的概率

```
k = 2, T = 240, t = 3, overlap_prob = 2.6%
k = 3, T = 240, t = 3, overlap_prob = 6.8%
k = 4, T = 240, t = 3, overlap_prob = 13.7%
k = 5, T = 240, t = 3, overlap_prob = 23.1%
k = 6, T = 240, t = 3, overlap_prob = 31.9%
k = 7, T = 240, t = 3, overlap_prob = 42.1%
k = 8, T = 240, t = 3, overlap_prob = 52.0%
k = 9, T = 240, t = 3, overlap_prob = 60.6%
k = 10, T = 240, t = 3, overlap_prob = 69.3%
k = 11, T = 240, t = 3, overlap_prob = 76.7%
k = 12, T = 240, t = 3, overlap_prob = 82.9%
k = 13, T = 240, t = 3, overlap_prob = 87.4%
k = 14, T = 240, t = 3, overlap_prob = 91.2%
k = 15, T = 240, t = 3, overlap_prob = 94.4%
k = 16, T = 240, t = 3, overlap_prob = 96.5%
k = 17, T = 240, t = 3, overlap_prob = 97.7%
k = 18, T = 240, t = 3, overlap_prob = 98.4%
k = 19, T = 240, t = 3, overlap_prob = 99.2%
k = 20, T = 240, t = 3, overlap_prob = 99.5%
k = 21, T = 240, t = 3, overlap_prob = 99.9%
k = 22, T = 240, t = 3, overlap_prob = 99.9%
k = 23, T = 240, t = 3, overlap_prob = 99.9%
k = 24, T = 240, t = 3, overlap_prob = 100.0%
k = 25, T = 240, t = 3, overlap_prob = 100.0%
```

利用递推公式验算$k=3,4$
```
from sympy import integrate, symbols

t = 3
T = 240
x = symbols('x')

P_B3 = 1-(2*t-T)**3 * (2*t+T) / (t-T) / T**3
P_B4 = 1-integrate(4*x**3 / T**4 * (3*t-x)**3 * (t+x) / (2*t-x) / (x-t)**3, (x, 3*t, T))

print("P_B3 = %.2f%%" % (P_B3*100))
print("P_B4 = %.2f%%" % (P_B4*100))
```

结果为

```
P_B3 = 3.79%
P_B4 = 5.18%
```

而蒙特卡洛模拟出来的结果是

```
k = 3, T = 240, t = 3, overlap_prob = 6.8%
k = 4, T = 240, t = 3, overlap_prob = 13.7%
```

也许哪里想错了？我再思考一下，和同学和教授再讨论一下吧
