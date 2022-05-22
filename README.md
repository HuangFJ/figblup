# FigBLUP
## 概述
品种选育在提高畜禽生产性能的诸多因素中作用最大，占40%左右。BLUP（最佳线性无偏预测）是评估动物遗传育种值的有效方法，自首次提出以来已经过去70多年。近年来，随着计算机的发展进步，其计算效率不断提升。目前BLUP已被许多国家采用，作为经济动物遗传评估的主要方法。

本项目主要目的是对BLUP进行求值，用Rust语言开发，力求计算高效，并顺便解答下面几个问题：
- 我在学习BLUP的理论知识，如何更好地去理解？
- 我用的BLUP不开源，它是否按照预设的参数和过程计算？
- 我想自己设计BLUP模型，有没有可以用来编写的基础原语？

## 基础
根据经典的遗传理论，个体的表型（phenotype）观测值是由环境因素和遗传因素共同决定的，具有以下简化模型：
$$\boldsymbol{u}=\boldsymbol{y} - \mu\boldsymbol{1} + \boldsymbol{e}$$
- $\boldsymbol{u}$是个体加性遗传值（即我们想要知道的真实育种值）且$\bar{u}=0$；
- $\boldsymbol{y}$是个体表型值，$\mu$是同一管理群的表型平均值（固定效应值）（已知）；
- $\boldsymbol{e}$是个体残差项且$\bar{e}=0$；
- $\boldsymbol{u}$和$\boldsymbol{e}$不相关$COV(u,e)=0$。

### 预测育种值
只考虑个体单一记录的情况下，预测育种值的一元回归方程为：
$$
\hat{\boldsymbol{u}}
=b(\boldsymbol{y}-\mu\boldsymbol{1})    \hspace8ex (1)
$$
- $\hat{\boldsymbol{u}}$即个体预测育种值；
- $\boldsymbol{y}-\mu\boldsymbol{1}$是零均值化（中心化）的个体表型值，$b$是其加权系数；
- 通过最小化$(\boldsymbol{u-\hat{u}})^2$求得$\boldsymbol{b}$：
$$
b
=\frac{COV(u,y)}{VAR(y)}  
=\frac{\sigma_u^2}{\sigma_y^2}
\Rightarrow h^2
$$

### 选择指数
如果考虑m组个体关系，此时预测育种值也叫选择指数，使用多元线性回归表示为：
$$
\hat{\boldsymbol{u}} 
=(Y_{n \times m} - \mu \boldsymbol{1})\boldsymbol{b}    \hspace8ex (2)
$$
- $Y_{n \times m} - \mu \boldsymbol{1}$是同一管理群相同表型的m组n个中心化观测数据；
- $\boldsymbol{b}=m \times 1$ 组的加权系数向量；
- 通过最小化$(\boldsymbol{u-\hat{u}})^2$得到$\boldsymbol{b}$的正规方程组：
$$
\boldsymbol{b}
=\frac{COV(u,y)}{VAR(y)} 
=\frac{(Y_{n \times m} - \mu \boldsymbol{1})^T\hat{u}}{(Y_{n \times m} - \mu \boldsymbol{1})^T(Y_{n \times m} - \mu \boldsymbol{1})}
\Rightarrow  \frac{G}{P}
$$
- $P$ 是每组个体的表型值的方差协方差矩阵，对角线元素为 $\sigma_y^2$；
- $G$ 是每组个体的表型值与预测育种值的方差协方差矩阵，对角线元素为 $\sigma_u^2=h^2\sigma_y^2$。


### 分子亲缘相关矩阵
个体间的加性遗传相关矩阵$A$，其主对角线元素$a_{ii}=1+F_i$，$F_i$是个体i的近交系数。$A$与加性遗传方差相乘$A\sigma_u^2$是个体育种值间的协方差，如果$u_i$是个体i的育种值，则$VAR(u_i)=a_{ii}\sigma_u^2$。

计算时，系谱中的个体从亲代到子代按1到n编号，按照以下原则计算：
- 如果个体i的双亲s和d已知，则：
    - $a_{ij}=a_{ji}=0.5(a_{js}+a_{jd});j=1,...,(i-1)$
    - $a_{ii}=1+0.5(a_{sd})$
- 如果仅有一个亲本s已知，并假设与配偶不相关：
    - $a_{ij}=a_{ji}=0.5(a_{js});j=1,...,(i-1)$
    - $a_{ii}=1$
- 如果双亲均未知，且双亲不相关：
    - $a_{ij}=a_{ji}=0;j=1,...,(i-1)$
    - $a_{ii}=1$


## 混合线性模型（BLUP）
综合选择指数和分子亲缘关系矩阵A，可以设计一个混合线性模型：
$$
\boldsymbol{y}=X\boldsymbol{\mu}+Z\boldsymbol{u}+\boldsymbol{e}
$$
- $\boldsymbol{y}=n \times 1$ 观测值向量，n是数据个数；
- $\boldsymbol{\mu}=p \times 1$ 固定效应向量，p是固定效应水平数；
- $\boldsymbol{u}=q \times 1$ 随机个体效应向量，q是随机效应水平数；
- $\boldsymbol{e}=n \times 1$ 随机残差效应向量；
- $X=n \times p$ 设计矩阵（incedence）；
- $Z=n \times q$ 设计矩阵（incedence）；
- $E(y)=X\mu,E(u)=E(e)=0$ ；
- $VAR(e)=I\sigma_e^2=R$ ；
- $VAR(u)=A\sigma_u^2=G$ ；
- $COV(u,e)=COV(e,u)=0$ ；
- $VAR(y)=V=ZGZ^T+R$ ；
- $COV(u,y)=GZ^T$ ；
- $COV(e,y)=R$ ;
- 用$y$的线性函数预测关于$\mu$和$u$的线性模型：
$$
k'\hat{\mu}+\hat{u}
=L'\boldsymbol{y}
$$
- 当$\hat{\boldsymbol{u}}=0$时：
$$
\hat{\boldsymbol{\mu}}
=(X^TV^{-1}X)X^TV^{-1}\boldsymbol{y}
$$
- 当$k'\hat{\mu}=0$时：
$$
\hat{\boldsymbol{u}}
=\boldsymbol{b} (\boldsymbol{y}-X\hat{\boldsymbol{\mu}})     \hspace8ex (3)
\\
$$

$$
\boldsymbol{b}
=\frac{COV(u,y)}{VAR(y)}
\Rightarrow  \frac{GZ^T}{V}
$$

### 选择指数代替个体的育种值：
- 用$\hat{\boldsymbol{\mu}}$代替选择指数中的$\mu$，BLUP与选择指数等价；
- 选择指数代替个体育种值：
$$
\hat{\boldsymbol{u}}=(I+\alpha A^{-1})^{-1} (\boldsymbol{y}-\hat{\boldsymbol{\mu}}\boldsymbol{1})
$$
- $Z=I;$
- $P=V=I \sigma_e^2 + A \sigma_u^2$
- $G=A \sigma_u^2$
- $\alpha=\sigma_e^2 / \sigma_u^2$ 或者 $\alpha=(1-h^2) / h^2$

### MME
综合以上所有内容：
$$
\begin{bmatrix}
X^TR^{-1}X & X^TR^{-1}Z \cr
Z^TR^{-1}X & Z^TR^{-1}Z+G^{-1}
\end{bmatrix}
\begin{bmatrix}
\hat{\mu} \cr
\hat{u}
\end{bmatrix}
=
\begin{bmatrix}
X^TR^{-1} y \cr
Z^TR^{-1} y
\end{bmatrix}
$$
因为$R$是单位矩阵，两边可以消去$R^{-1}$，这就是我们需要的全部：
$$
\begin{bmatrix}
X^TX & X^TZ \cr
Z^TX & Z^TZ+\alpha A^{-1}
\end{bmatrix}
\begin{bmatrix}
\hat{\mu} \cr
\hat{u}
\end{bmatrix}
=
\begin{bmatrix}
X^T y \cr
Z^T y
\end{bmatrix}
\hspace8ex (4)
$$
