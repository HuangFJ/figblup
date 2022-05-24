# FigBLUP
## 概述
线性混合模型既不需要假设固定效应为已知，也不需要假设观测数据的独立性和方差齐性，适合处理平衡/不平衡的纵向数据（Panel Data）、重复测量数据、平衡/不平衡试验设计数据、多环境试验数据、区组数据以及空间相关数据。而一元和多元线性模型则需要假设固定效应为已知，而且前者方差协方差矩阵只是一个标量，后者方差协方差矩阵是无结构的，线性混合模型可以根据数据本身的结构特点，较为灵活地选择其方差协方差矩阵结构。因此，线性混合模型广泛地应用在生物、医学、经济、金融、环境科学及抽样调查等领域。

本文介绍线性混合模型在动物中的应用，以 $A$ 作为方差协方差矩阵结构，对其随机效应（遗传育种值）进行 BLUP（最佳线性无偏预测）分析。这个工作非常有意义，因为在提高畜禽生产性能的诸多因素中，遗传评估占 40% 左右，作用是最大的。

本文源码用 Rust 语言编写，力求计算高效，并顺便解答下面几个问题：
- 一元/多元线性模型和线性混合模型各有什么特点？
- 我用的 BLUP 不开源，它是怎么实现的？
- 我想在其他场景中应用线性混合模型及 BLUP，有没有可以用来参考的案例？

## 基础
根据经典的遗传理论，个体的表型（phenotype）观测值是由环境因素（固定效应）和遗传因素（随机效应）共同决定的，具有以下简化模型：
$$\boldsymbol{y}=\mu\boldsymbol{1}+\boldsymbol{u}+\boldsymbol{e}$$
> $\boldsymbol{e}\sim\mathcal{N}(0,\sigma_e^2)$ 是个体残差项；<br>
> $\boldsymbol{u}\sim\mathcal{N}(0,\sigma_u^2)$ 是个体加性遗传值（即我们想要知道的真实育种值）；<br>
> $\boldsymbol{y}\sim\mathcal{N}(\mu,\sigma_y^2)$ 是个体表型值，$\mu$ 是已知同一管理群的表型平均值；<br>
> $COV(u,e)=COV(e,u)=0$ 即 $\boldsymbol{u}$ 和 $\boldsymbol{e}$ 不相关。<br>

### 一元线性模型
只考虑个体单一记录的情况下，个体预测育种值 $\hat{\boldsymbol{u}}$ 的一元回归方程为：
$$\hat{\boldsymbol{u}}=b(\boldsymbol{y}-\mu\boldsymbol{1})    \hspace8ex (1)$$

通过最小化 $(\boldsymbol{u-\hat{u}})^2$ 求得 $b$：
$$b=\frac{COV(u,y)}{VAR(y)}=\frac{\sigma_u^2}{\sigma_y^2} \Rightarrow h^2$$

### 多元线性模型
如果考虑 $m$ 组记录，此时个体预测育种值也叫选择指数，使用多元线性回归方程为：
$$\hat{\boldsymbol{u}}=(Y_{n \times m}-\boldsymbol{\mu}^T\otimes\boldsymbol{1})\boldsymbol{b}    \hspace8ex (2)$$
> $\boldsymbol{\mu}=m \times 1$ 是已知各组的平均值向量；<br>
> $\boldsymbol{b}=m \times 1$ 组的加权系数向量。<br>

通过最小化 $(\boldsymbol{u-\hat{u}})^2$ 得到 $\boldsymbol{b}$ 的正规方程组解：
$$\boldsymbol{b}=\frac{COV(u,y)}{VAR(y)}=\frac{(Y_{n \times m}-\boldsymbol{\mu}^T\otimes\boldsymbol{1})^T\hat{u}}{(Y_{n \times m}-\boldsymbol{\mu}^T\otimes\boldsymbol{1})^T(Y_{n \times m}-\boldsymbol{\mu}^T\otimes\boldsymbol{1})} \Rightarrow \frac{G}{P}$$
> $P=m \times m$ 是观测值的方差协方差矩阵，对角线元素为表型方差 $\sigma_y^2$ ；<br>
> $G=m \times 1$ 是观测值与预测育种值间的协方差矩阵，对角线元素为遗传方差 $\sigma_u^2=h^2\sigma_y^2$ 。<br>

### 分子亲缘相关矩阵
个体间的加性遗传相关矩阵 $A$，其主对角线元素 $a_{ii}=1+F_i$，$F_i$ 是个体 $i$ 的近交系数。$A$ 与加性遗传方差相乘 $A\sigma_u^2$ 是个体育种值间的协方差，如果 $u_i$ 是个体 $i$ 的育种值，则 $VAR(u_i)=a_{ii}\sigma_u^2$ 。

计算时，系谱中的个体从亲代到子代按 $1$ 到 $n$ 编号，按照以下原则计算：
> 1、如果个体 $i$ 的双亲 $s$ 和 $d$ 已知，则：
>> $a_{ij}=a_{ji}=0.5(a_{js}+a_{jd});j=1,...,(i-1)$ <br>
>> $a_{ii}=1+0.5(a_{sd})$
>
> 2、如果仅有一个亲本 $s$ 已知，并假设与配偶不相关：
>> $a_{ij}=a_{ji}=0.5(a_{js});j=1,...,(i-1)$ <br>
>> $a_{ii}=1$
>
> 3、如果双亲均未知，且双亲不相关：
>> $a_{ij}=a_{ji}=0;j=1,...,(i-1)$ <br>
>> $a_{ii}=1$

## 线性混合模型
$$\boldsymbol{y}=X\boldsymbol{\mu}+Z\boldsymbol{u}+\boldsymbol{e}$$
> $\boldsymbol{e}\sim\mathcal{N}(0,R=I\sigma_e^2)$ 是一个 $n \times 1$ 随机残差效应向量，$n$ 是数据个数；<br>
> $\boldsymbol{u}\sim\mathcal{N}(0,G=A\sigma_u^2)$ 是一个 $q \times 1$ 的随机个体效应向量，$q$ 是随机效应水平数；<br>
> $\boldsymbol{y}\sim\mathcal{N}(X\mu,V)$ 是一个 $n \times 1$ 的观测值向量，$n$ 是数据个数；
>> $X=n \times p$ 设计矩阵（incedence）；<br>
>> $\boldsymbol{\mu}=p \times 1$ 固定效应向量，$p$ 是固定效应水平数；<br>
>> $Z=n \times q$ 设计矩阵（incedence）；<br>
>> $V=ZGZ^T+R$ ；<br>
>
> $COV(u,e)=COV(e,u)=0$ 即 $\boldsymbol{u}$ 和 $\boldsymbol{e}$ 不相关；<br>
> $COV(u,y)=GZ^T$ ；<br>
> $COV(e,y)=R$ 。

用 $y$ 的线性函数预测关于 $\mu$ 和 $u$ 的线性模型：
$$k'\hat{\mu}+\hat{u}=L'\boldsymbol{y}$$

当 $\hat{\boldsymbol{u}}=0$ 时，BLUE：
$$\hat{\boldsymbol{\mu}}=(X^TV^{-1}X)^{-1}X^TV^{-1}\boldsymbol{y}$$

当 $k'\hat{\mu}=0$ 时，BLUP：
$$\hat{\boldsymbol{u}}=\boldsymbol{b} (\boldsymbol{y}-X\hat{\boldsymbol{\mu}})     \hspace8ex (3)$$

$$\boldsymbol{b}=\frac{COV(u,y)}{VAR(y)} \Rightarrow  \frac{GZ^T}{V}$$

### 选择指数代替个体的育种值：
用 $\hat{\boldsymbol{\mu}}$ 代替选择指数中的 $\boldsymbol{\mu}$，则 BLUP 与选择指数等价。特定情况下，选择指数可以代替个体育种值：
$$\hat{\boldsymbol{u}}=(I+\alpha A^{-1})^{-1} (\boldsymbol{y}-X\hat{\boldsymbol{\mu}})$$
> $Z=I$ <br>
> $G=A \sigma_u^2$ <br>
> $P=V=A \sigma_u^2 + I \sigma_e^2$ <br>
> $\alpha=\sigma_e^2 / \sigma_u^2$ 或者 $\alpha=(1-h^2) / h^2$ <br>

### MME
如果假设一些前提条件，则可以用 MME 计算 BLUE 和 BLUP 值：
> $\boldsymbol{e}, \boldsymbol{u}, \boldsymbol{y}$ 分布为多元正态分布，$\sigma_y^2=\sigma_u^2 + \sigma_e^2$。遗传方差变化是因为选择（比如选育近交）导致的，这个变化用 $A$ 矩阵体现；<br>
> 基础群的 $\sigma_u^2$ 和 $\sigma_e^2$ 是已知的。但实际上未知，它们可以用约束最大似然 REML 进行估计。<br>

$$\begin{bmatrix} X^TR^{-1}X & X^TR^{-1}Z \cr Z^TR^{-1}X & Z^TR^{-1}Z+G^{-1} \end{bmatrix}\begin{bmatrix}\hat{\mu} \cr \hat{u} \end{bmatrix}=\begin{bmatrix}X^TR^{-1} y \cr Z^TR^{-1} y \end{bmatrix}$$

因为 $R$ 是单位矩阵，两边可以消去 $R^{-1}$：
$$\begin{bmatrix} X^TX & X^TZ \cr Z^TX & Z^TZ+\alpha A^{-1} \end{bmatrix} \begin{bmatrix}\hat{\mu} \cr \hat{u} \end{bmatrix}=\begin{bmatrix} X^T y \cr Z^T y \end{bmatrix} \hspace8ex (4)$$

## 约束最大似然
